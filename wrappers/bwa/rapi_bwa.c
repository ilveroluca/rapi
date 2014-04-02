/*
 * rapi_bwa.c
 */

#include "aligner.h"


aln_kv* aln_kv_insert_char(const char* key, char value,           const aln_kv* tail);
aln_kv* aln_kv_insert_str( const char* key, const char* value,    const aln_kv* tail);
aln_kv* aln_kv_insert_long(const char* key, long value,           const aln_kv* tail);
aln_kv* aln_kv_insert_dbl( const char* key, double value,         const aln_kv* tail);
aln_kv* aln_kv_insert_ba(  const char* key, const uint8_t* value, const aln_kv* tail);
aln_kv* aln_kv_insert_la(  const char* key, const long* value,    const aln_kv* tail);

aln_kv_free_list(aln_kv* list);



/* Init Library */
int aln_init()
{
  // DO

  return ALN_NO_ERROR;
}

/* Init Library Options */
int aln_init_ops( aln_opts * my_opts )
{
 
}

/* Load Reference */
const char * aln_version()
{
 
}

/* Load Reference */
int aln_load_ref( const char * reference, aln_ref * ref_struct )
{
  ref_struct->handle=bwa_idx_load(reference, BWA_IDX_ALL);
  if ( ref_struct->handle == 0 )
  {
    return ALN_REFERENCE_ERROR;
  }
 
  ref_struct->path=reference;
 
  /* Fill in Contig Information */
  ref_struct->n_contigs = reference->handle->bns->n_seqs;	/* Handle contains bnt_seq_t * bns holding contig information */
  ref_struct->contigs=calloc( ref_struct->n_contigs, sizeof(aln_contig) );
 
  if ( ref_struct->contigs == 0 )
  {
    return ALN_MEMORY_ERROR;
  }
 
  for ( int contig_loop=0, contig_loop < ref_struct->n; contig_loop++ )
  {
    /* Get Name */
    /* Each bns holds an array of bntann1_t *anns -> holding 	
     	int64_t offset;
	int32_t len;
	int32_t n_ambs;
	uint32_t gi;
	char *name, *anno;
    */
    ref_struct->contigs[ contig_loop ].name =  reference->handle->bns->anns[ contig_loop ].name;
    ref_struct->contigs[ contig_loop ].len =  reference->handle->bns->anns[ contig_loop ].len;

  }

  return ALN_NO_ERROR;
}

/* Free Reference */
int aln_free_ref( aln_ref * ref_struct )
{
  for ( int contig_loop=0, contig_loop < ref_struct->n; contig_loop++ )
  {
    free (  ref_struct->contigs );
  }
}

/* Allocate reads */
int aln_alloc_reads( aln_batch * batch, int n_reads_fragment, int n_fragments );
{
    batch=calloc( n_reads_fragment * n_fragments, sizeof(aln_batch) );
    return ( batch == 0 ? ALN_MEMORY_ERROR : ALN_NO_ERROR );
}

int aln_set_read(const aln_batch* batch, int n_frag, int n_read, const char* name, const char* seq, const char* qual, int q_offset) {
	int error_code = 0;

	if (n_frag >= batch->n_frags || n_read >= batch->n_reads_frag)
		return ALN_PARAM_ERROR;

	aln_read = batch->reads[n_frag * batch->n_reads_frags + n_read];
	const int name_len = strlen(name);
	const int seq_len = strlen(seq);
	seq->length = seq->clip_len = seq_len;

	seq->seq = (ubyte_t*)malloc(seq_len);
	seq->name = (char*)malloc(name_len+1); // +1 for \0
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
	// j
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
