/*
 * rapi_bwa.c
 */

#include <aligner.h>
#include <bwamem.h>

#include <string.h>
#include <stdio.h>
#include <stdlib.h>

char* bwa_pg = "rapi";

#if 0
aln_kv* aln_kv_insert_char(const char* key, char value,           const aln_kv* tail);
aln_kv* aln_kv_insert_str( const char* key, const char* value,    const aln_kv* tail);
aln_kv* aln_kv_insert_long(const char* key, long value,           const aln_kv* tail);
aln_kv* aln_kv_insert_dbl( const char* key, double value,         const aln_kv* tail);
aln_kv* aln_kv_insert_ba(  const char* key, const uint8_t* value, const aln_kv* tail);
aln_kv* aln_kv_insert_la(  const char* key, const long* value,    const aln_kv* tail);

aln_kv_free_list(aln_kv* list);
#endif

#if 1 // strdup is not defined if we compile with c99
char* strdup(const char* str)
{
	if (NULL == str)
		return NULL;

	size_t len = strlen(str);
	char* new_str = malloc(len + 1);
	if (NULL == new_str)
		return NULL;
	return strcpy(new_str, str);
}
#endif

/* Init Library */
int aln_init(const aln_opts* opts)
{
	mem_opt_t*const bwa_opt = (mem_opt_t*) opts->_private;

	bwa_fill_scmat(bwa_opt->a, bwa_opt->b, bwa_opt->mat);

	return ALN_NO_ERROR;
}

/* Init Library Options */
int aln_init_opts( aln_opts * my_opts )
{
	// create a BWA opt structure
	mem_opt_t*const bwa_opt = mem_opt_init();

	// set the default values as in BWA fastmap.c
	bwa_opt->a = bwa_opt->b = bwa_opt->o_del = bwa_opt->e_del = -1;
	bwa_opt->o_ins = bwa_opt->e_ins = bwa_opt->pen_unpaired = -1;
	bwa_opt->pen_clip5 = bwa_opt->pen_clip3 = bwa_opt->zdrop = bwa_opt->T = -1;

	my_opts->_private = bwa_opt; // hand the bwa structure onto the external one
	my_opts->ignore_unsupported = 1;
	my_opts->mapq_min     = 0;
	my_opts->trim_quality = 0;
	my_opts->isize_min    = 0;
	my_opts->isize_max    = bwa_opt->max_ins;
	my_opts->n_parameters = 0;
	my_opts->parameters   = NULL;

	return ALN_NO_ERROR;
}

int aln_free_opts( aln_opts * my_opts )
{
	free(my_opts->_private);
	return ALN_NO_ERROR;
}

/* Load Reference */
const char * aln_version()
{
	return "my version!";
}

/* Load Reference */
int aln_load_ref( const char * reference_path, aln_ref * ref_struct )
{
	if ( NULL == ref_struct || NULL == reference_path )
		return ALN_PARAM_ERROR;

	const bwaidx_t*const bwa_idx = bwa_idx_load(reference_path, BWA_IDX_ALL);
	if ( NULL == bwa_idx )
		return ALN_REFERENCE_ERROR;

	// allocate memory
	ref_struct->path = strdup(reference_path);
	ref_struct->contigs = calloc( bwa_idx->bns->n_seqs, sizeof(aln_contig) );
	if ( NULL == ref_struct->path || NULL == ref_struct->contigs )
	{
		// if either allocations we free everything and return an error
		bwa_idx_destroy((bwaidx_t*)bwa_idx);
		free(ref_struct->path);
		free(ref_struct->contigs);
		memset(ref_struct, 0, sizeof(*ref_struct));
		return ALN_MEMORY_ERROR;
	}
 
	/* Fill in Contig Information */
	ref_struct->n_contigs = bwa_idx->bns->n_seqs; /* Handle contains bnt_seq_t * bns holding contig information */
	for ( int i = 0; i < ref_struct->n_contigs; ++i )
	{
		ref_struct->contigs[ i ].name = bwa_idx->bns->anns[ i ].name;
		ref_struct->contigs[ i ].len = bwa_idx->bns->anns[ i ].len;
	}
	ref_struct->_private = (bwaidx_t*)bwa_idx; // cast to discard const

	return ALN_NO_ERROR;
}

/* Free Reference */
int aln_free_ref( aln_ref * ref )
{
	// free bwa's part
	bwa_idx_destroy(ref->_private);

	// then free the rest of the structure
	free(ref->path);

	if (ref->contigs) {
		for ( int i = 0; i < ref->n_contigs; ++i )
		{
			aln_contig* c = &ref->contigs[i];
			free(c->name);
			free(c->assembly_identifier);
			free(c->species);
			free(c->uri);
		}
		free (ref->contigs);
	}
	memset(ref, 0, sizeof(*ref));
	return ALN_NO_ERROR;
}

/* Allocate reads */
int aln_alloc_reads( aln_batch * batch, int n_reads_fragment, int n_fragments )
{
	if (n_fragments < 0 || n_reads_fragment < 0)
		return ALN_PARAM_ERROR;

	batch->reads = calloc( n_reads_fragment * n_fragments, sizeof(aln_read) );
	if (NULL == batch->reads)
		return ALN_MEMORY_ERROR;
	batch->n_frags = n_fragments;
	batch->n_reads_frag = n_reads_fragment;
	return ALN_NO_ERROR;
}

int aln_free_reads( aln_batch * batch )
{
	for (int f = 0; f < batch->n_frags; ++f) {
		for (int r = 0; r < batch->n_reads_frag; ++r) {
			aln_read* read = aln_get_read(batch, f, r);
			free(read->id);
			free(read->seq);
			free(read->qual);
		}
	}

	free(batch->reads);
	memset(batch, 0, sizeof(*batch));

	return ALN_NO_ERROR;
}

int aln_set_read(aln_batch* batch,
			int n_frag, int n_read,
			const char* name, const char* seq, const char* qual,
			int q_offset) {
	int error_code = 0;

	if (n_frag >= batch->n_frags || n_read >= batch->n_reads_frag)
		return ALN_PARAM_ERROR;

	aln_read* read = aln_get_read(batch, n_frag, n_read);
	const int name_len = strlen(name);
	const int seq_len = strlen(seq);
	read->length = read->clip_len = seq_len;

	read->seq = (char*)malloc(seq_len);
	read->id = (char*)malloc(name_len + 1); // +1 for \0
	if (qual)
		read->qual = (ubyte_t*)malloc(seq_len);
	else
		read->qual = NULL;

	if (read->seq == NULL || read->id == NULL || (read->qual == NULL && qual != NULL))
	{
		fprintf(stderr, "Unable to allocate memory for sequence\n");
		error_code = ALN_MEMORY_ERROR;
		goto error;
	}
	// else

	// encode sequence and set quality
	/** XXX: should we recode the bases here? */
	memcpy(read->seq, seq, seq_len);

	if (read->qual)
	{
		for (int i = 0; i < seq_len; ++i)
		{
			read->qual[i] = (int)qual[i] - q_offset + 33; // 33 is the Sanger offset.  BWA expects it this way.
			if (read->qual[i] > 127)
			{ // qual is unsigned, by Sanger base qualities have an allowed range of [0,94], and 94+33=127
				fprintf(stderr, "Invalid base quality score %d\n", read->qual[i]);
				error_code = ALN_PARAM_ERROR;
				goto error;
			}
		}
	}

	// and finally the name
	strcpy(read->id, name);
	// trim /[12]$
	int t = name_len;
	if (t > 2 && read->id[t-2] == '/' && (read->id[t-1] == '1' || read->id[t-1] == '2'))
		read->id[t-2] = '\0';
	return 0;

error:
	// In case of error, free any allocated memory and return -1
	free(read->seq); free(read->id); free(read->qual);
	return error_code;
}

int aln_align_reads( const aln_ref* ref,  aln_batch * batch, const aln_opts * config )
{
	/*
	mem_opt_t*const bwa_opt = (mem_opt_t*) config->_private;
	
	const bwaidx_t*const bwa_idx = ref->_private;

	mem_alnreg_v bwa_alignment;
	*/

	// traslate our read structure into BWA reads

	// run the alignment
	//bwa_alignment = mem_align1(bwa_opt, bwa_idx->bwt, bwa_idx->bns, bwa_idx->pac, ks->seq.l, ks->seq.s); // get all the hits
	// translate BWA's results back into our structure.
	return ALN_NO_ERROR;
}

