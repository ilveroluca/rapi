/*
 * rapi_bwa.c
 */

#include "../../aligner.h"
#include <bwamem.h>

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/**
 * The 'bwa_pg' string is statically allocated in some BWA file that we're not
 * linking, so we need to define it here.  I think it's the string that the
 * program uses to identify itself in the SAM header's @PG tag.
 */
char* bwa_pg = "rapi";

/******** Utility functions *******/
void aln_print_read(FILE* out, const aln_read* read)
{
	fprintf(out, "read id: %s\n", read->id);
	fprintf(out, "read length: %d\n", read->length);
	fprintf(out, "read seq: %.*s\n", read->length, read->seq);
	fprintf(out, "read qual: not implemented\n");
}
/**********************************/

/**
 * Definition of the aligner state structure.
 */
struct aln_aligner_state {
	int64_t n_reads_processed;
	// paired-end stats
	mem_pestat_t pes[4];
	mem_pestat_t pes0;
};

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
	// no op
	return ALN_NO_ERROR;
}

/* Init Library Options */
int aln_init_opts( aln_opts * my_opts )
{
	// create a BWA opt structure
	mem_opt_t*const bwa_opt = mem_opt_init();

	// Default values copied from bwamem.c in 0.7.8
	bwa_opt->flag = 0;
	bwa_opt->a = 1; bwa_opt->b = 4;
	bwa_opt->o_del = bwa_opt->o_ins = 6;
	bwa_opt->e_del = bwa_opt->e_ins = 1;
	bwa_opt->w = 100;
	bwa_opt->T = 30;
	bwa_opt->zdrop = 100;
	bwa_opt->pen_unpaired = 17;
	bwa_opt->pen_clip5 = bwa_opt->pen_clip3 = 5;
	bwa_opt->min_seed_len = 19;
	bwa_opt->split_width = 10;
	bwa_opt->max_occ = 10000;
	bwa_opt->max_chain_gap = 10000;
	bwa_opt->max_ins = 10000;
	bwa_opt->mask_level = 0.50;
	bwa_opt->chain_drop_ratio = 0.50;
	bwa_opt->split_factor = 1.5;
	bwa_opt->chunk_size = 10000000;
	bwa_opt->n_threads = 1;
	bwa_opt->max_matesw = 100;
	bwa_opt->mask_level_redun = 0.95;
	bwa_opt->mapQ_coef_len = 50; bwa_opt->mapQ_coef_fac = log(bwa_opt->mapQ_coef_len);
	bwa_fill_scmat(bwa_opt->a, bwa_opt->b, bwa_opt->mat);

	my_opts->_private = bwa_opt; // hand the bwa structure onto the external one
	my_opts->ignore_unsupported = 1;
	my_opts->mapq_min     = 0;
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


typedef struct {
	unsigned long n_bases;
	int n_reads;
	int n_reads_per_frag;
	bseq1_t* seqs;
} bwa_batch;

static void _free_bwa_batch_contents(bwa_batch* batch)
{
	for (int i = 0; i < batch->n_reads; ++i)
		free(batch->seqs[i].seq);

	free(batch->seqs);
	memset(batch, 0, sizeof(bwa_batch));
}


static int _batch_to_bwa_seq(const aln_batch* batch, const aln_opts* opts, bwa_batch* bwa_seqs)
{
	bwa_seqs->n_bases = 0;
	bwa_seqs->n_reads = 0;
	bwa_seqs->n_reads_per_frag = batch->n_reads_frag;
	fprintf(stderr, "Need to allocate %d elements of size %ld\n", batch->n_frags * batch->n_reads_frag, sizeof(bseq1_t));

	bwa_seqs->seqs = calloc(batch->n_frags * batch->n_reads_frag, sizeof(bseq1_t));
	if (NULL == bwa_seqs->seqs) {
		fprintf(stderr, "Allocation failed!\n");
		return ALN_MEMORY_ERROR;
	}

	for (int f = 0; f < batch->n_frags; ++f)
	{
		for (int r = 0; r < batch->n_reads_frag; ++r)
		{
			const aln_read* rapi_read = aln_get_read(batch, f, r);
			bseq1_t* bwa_read = bwa_seqs->seqs + bwa_seqs->n_reads;

			int buf_size = rapi_read->length;
			// We place the seq and qual arrays in a single allocated buffer.
			// It's important to remember this when freeing the structure. On the
			// other hand, the other strings are not duplicated since BWA doesn't
			// modify them.
			if (rapi_read->qual)
				buf_size *= 2;
			bwa_read->seq = malloc(buf_size);
			if (NULL == bwa_read->seq)
				goto failed_allocation;
			memcpy(bwa_read->seq, rapi_read->seq, rapi_read->length);
			if (rapi_read->qual) {
				bwa_read->qual = bwa_read->seq + rapi_read->length;
				memcpy(bwa_read->qual, rapi_read->qual, rapi_read->length);
			}
			else
				bwa_read->qual = NULL;

			bwa_read->name = rapi_read->id;
			bwa_read->l_seq = rapi_read->length;
			bwa_read->comment = NULL;
			bwa_read->sam = NULL;

			// As we build the objects n_reads keeps the actual number
			// constructed thus far and thus can be used to free the allocated
			// structres in case of error.
			bwa_seqs->n_reads += 1;
			bwa_seqs->n_bases += rapi_read->length;
		}
	}
	return ALN_NO_ERROR;

failed_allocation:
	fprintf(stderr, "Failed to allocate while constructing sequences! Freeing and returning\n");
	_free_bwa_batch_contents(bwa_seqs);
	return ALN_MEMORY_ERROR;
}

/*
 * Many of the default option values need to be adjusted if the matching score
 * (opt->a) is changed.  This function (from the BWA code) does that.
 *
 * \param opt The structure containing actual values and that will be "fixed"
 * \param override set the members to 1 if the value has been overridden and thus should be kept;
 *                 else the corresponding value in opts is assumed to be at default and will be adjusted.
 *
 * So, if you change 'a' set the 'override' accordingly and call this function.  E.g.,
 *
 * <pre>
 *   aln_opts* opt = aln_init_opts();
 *   mem_opt_t* bwa_opts = (mem_opt_t*) opt->_private;
 *   bwa_opts->a = 2;
 *   bwa_opts->b = 5;
 *   mem_opt_t override;
 *   override.a = override.b = 1;
 *   adjust_bwa_opts(bwa_opts, &override);
 * </pre>
 */
void adjust_bwa_opts(mem_opt_t* opt, const mem_opt_t* override)
{
	if (override->a == 1) { // matching score is changed
		if (override->b != 1) opt->b *= opt->a;
		if (override->T != 1) opt->T *= opt->a;
		if (override->o_del != 1) opt->o_del *= opt->a;
		if (override->e_del != 1) opt->e_del *= opt->a;
		if (override->o_ins != 1) opt->o_ins *= opt->a;
		if (override->e_ins != 1) opt->e_ins *= opt->a;
		if (override->zdrop != 1) opt->zdrop *= opt->a;
		if (override->pen_clip5 != 1) opt->pen_clip5 *= opt->a;
		if (override->pen_clip3 != 1) opt->pen_clip3 *= opt->a;
		if (override->pen_unpaired != 1) opt->pen_unpaired *= opt->a;
		bwa_fill_scmat(opt->a, opt->b, opt->mat);
	}
}

static int _convert_opts(const aln_opts* opts, mem_opt_t* bwa_opts)
{
	bwa_opts->T = opts->mapq_min;
	bwa_opts->max_ins = opts->isize_max;

	// TODO: other options provided through 'parameters' field
	return ALN_NO_ERROR;
}

int aln_init_aligner_state(const aln_opts* opts, struct aln_aligner_state** ret_state)
{
	// allocate and zero the structure
	aln_aligner_state* state = *ret_state = calloc(1, sizeof(aln_aligner_state));
	if (NULL == state)
		return ALN_MEMORY_ERROR;
	return ALN_NO_ERROR;
}

int aln_free_aligner_state(aln_aligner_state* state)
{
	free(state);
	return ALN_NO_ERROR;
}

int aln_align_reads( const aln_ref* ref,  aln_batch * batch, const aln_opts * config, aln_aligner_state* state )
{
	int error = 0;

	if (batch->n_reads_frag > 2)
		return ALN_OP_NOT_SUPPORTED_ERROR;

	if (batch->n_reads_frag <= 0)
		return ALN_PARAM_ERROR;

	// "extract" BWA-specific structures
	mem_opt_t*const bwa_opt = (mem_opt_t*) config->_private;
	const bwaidx_t*const bwa_idx = ref->_private;

	if (batch->n_reads_frag == 2) // paired-end
		bwa_opt->flag |= MEM_F_PE;

	if ((error = _convert_opts(config, bwa_opt)))
		return error;
	fprintf(stderr, "opts converted\n");

	// traslate our read structure into BWA reads
	bwa_batch bwa_seqs;
	if ((error = _batch_to_bwa_seq(batch, config, &bwa_seqs)))
		return error;
	fprintf(stderr, "converted reads to BWA structures.\nGoing to process.\n");

	// run the alignment
	mem_process_seqs(bwa_opt, bwa_idx->bwt, bwa_idx->bns, bwa_idx->pac, state->n_reads_processed, bwa_seqs.n_reads, bwa_seqs.seqs, &state->pes0);
	state->n_reads_processed += bwa_seqs.n_reads;
	fprintf(stderr, "processed %ld reads\n", state->n_reads_processed);

	//mem_alnreg_v bwa_alignment;
	//bwa_alignment = mem_align1(bwa_opt, bwa_idx->bwt, bwa_idx->bns, bwa_idx->pac, ks->seq.l, ks->seq.s); // get all the hits
	// translate BWA's results back into our structure.

	_free_bwa_batch_contents(&bwa_seqs);

	return ALN_NO_ERROR;
}

/******* Generic functions ***********
 * These are generally usable. We should put them in a generic rapi.c file.
 *************************************/

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
	read->length = seq_len;

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

