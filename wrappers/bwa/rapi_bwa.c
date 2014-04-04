/*
 * rapi_bwa.c
 */

#include "aligner.h"
#include "bwamem.h"
#include <string.h>


aln_kv* aln_kv_insert_char(const char* key, char value,           const aln_kv* tail);
aln_kv* aln_kv_insert_str( const char* key, const char* value,    const aln_kv* tail);
aln_kv* aln_kv_insert_long(const char* key, long value,           const aln_kv* tail);
aln_kv* aln_kv_insert_dbl( const char* key, double value,         const aln_kv* tail);
aln_kv* aln_kv_insert_ba(  const char* key, const uint8_t* value, const aln_kv* tail);
aln_kv* aln_kv_insert_la(  const char* key, const long* value,    const aln_kv* tail);

aln_kv_free_list(aln_kv* list);


#if 0
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
int aln_init()
{
	// no op
  return ALN_NO_ERROR;
}

/* Init Library Options */
int aln_init_ops( aln_opts * my_opts )
{
	// create a BWA opt structure and link it to the external one
	my_opts->_private = mem_opt_init();

	my_opts->ignore_unsupported = TRUE;
	my_opts->mapq_min     = 0;
	my_opts->trim_quality = 0;
	my_opts->isize_min    = 0;
	my_opts->isize_max    = my_opts->_private.max_ins;
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

  ref_struct->_private = bwa_idx_load(reference, BWA_IDX_ALL);

  if ( NULL == ref_struct->_private )
    return ALN_REFERENCE_ERROR;
 
  ref_struct->path = strdup(reference_path);
 
  /* Fill in Contig Information */
  ref_struct->n_contigs = reference->handle->bns->n_seqs;	/* Handle contains bnt_seq_t * bns holding contig information */
  ref_struct->contigs = calloc( ref_struct->n_contigs, sizeof(aln_contig) );
 
  if ( NULL == ref_struct->contigs )
  {
    return ALN_MEMORY_ERROR;
  }
 
  for ( int contig_loop = 0, contig_loop < ref_struct->n_contigs; contig_loop++ )
  {
    /* Get Name */
    /* Each bns holds an array of bntann1_t *anns -> holding 	
     	int64_t offset;
	int32_t len;
	int32_t n_ambs;
	uint32_t gi;
	char *name, *anno;
    */
    ref_struct->contigs[ contig_loop ].name = reference->handle->bns->anns[ contig_loop ].name;
    ref_struct->contigs[ contig_loop ].len =  reference->handle->bns->anns[ contig_loop ].len;
  }

  return ALN_NO_ERROR;
}

/* Free Reference */
int aln_free_ref( aln_ref * ref_struct )
{
	free(ref_struct->path);

  for ( int contig_loop = 0, contig_loop < ref_struct->n_contigs; contig_loop++ )
  {
    free (  ref_struct->contigs );
  }
}

/* Allocate reads */
int aln_alloc_reads( aln_batch * batch, int n_reads_fragment, int n_fragments );
{
	if (n_fragments < 0 || n_reads_fragment < 0)
		return ALN_PARAM_ERROR;

	batch = calloc( n_reads_fragment * n_fragments, sizeof(aln_batch) );
	batch->n_fragments = n_fragments;
	batch->n_reads_fragment = n_reads_fragment;
	return ( batch == 0 ? ALN_MEMORY_ERROR : ALN_NO_ERROR );
}

int aln_free_reads( aln_batch * batch )
{
	for (int f = 0; f < batch->n_frags; ++f) {
		for (int r = 0; r < batch->n_reads_frags; ++r) {
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

int aln_set_read(const aln_batch* batch,
	    int n_frag, int n_read,
	    const char* name, const char* seq, const char* qual,
	    int q_offset) {
	int error_code = 0;

	if (n_frag >= batch->n_frags || n_read >= batch->n_reads_frag)
		return ALN_PARAM_ERROR;

	aln_read* seq = aln_get_read(batch, n_frag, n_read);
	const int name_len = strlen(name);
	const int seq_len = strlen(seq);
	seq->length = seq->clip_len = seq_len;

	seq->seq = (char*)malloc(seq_len);
	seq->name = (char*)malloc(name_len + 1); // +1 for \0
	if (qual)
		seq->qual = (ubyte_t*)malloc(seq_len);
	else
		seq->qual = NULL;

	if (seq->seq  == NULL || seq->name ==  NULL || (seq->qual ==  NULL && qual != NULL))
	{
		fprintf(stderr, "Unable to allocate memory for sequence\n");
		error_code = ALN_MEMORY_ERROR;
		goto error;
	}
	// else

	// encode sequence and set quality
	/** XXX: should we recode the bases here? */
	memcpy(seq->seq, seq, seq_len);

	if (seq->qual)
	{
		for (int i = 0; i < seq_len; ++i)
		{
			seq->qual[i] = (int)qual[i] - q_offset + 33; // 33 is the Sanger offset.  BWA expects it this way.
			if (seq->qual[i] > 127)
			{ // qual is unsigned, by Sanger base qualities have an allowed range of [0,94], and 94+33=127
				fprintf(stderr, "Invalid base quality score %d\n", seq->qual[i]);
				error_code = ALN_PARAM_ERROR;
				goto error;
			}
		}
	}

	// and finally the name
	strcpy(seq->id, name);
	// trim /[12]$
	int t = name_len;
	if (t > 2 && seq->id[t-2] == '/' && (seq->id[t-1] == '1' || seq->id[t-1] == '2'))
		seq->name[t-2] = '\0';
	return 0;

error:
	// In case of error, free any allocated memory and return -1
	free(seq->seq); free(seq->id); free(seq->qual);
	return error_code;
}
