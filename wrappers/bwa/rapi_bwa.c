/*
 * rapi_bwa.c 
 */

#include aligner.h
#include ... bwa headers ...




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
    return ( batch == 0 ? ALN_MEMORY_ERROR, ALN_NO_ERROR );
}

/* Set read */
int aln_set_read(aln_batch * batch, int n_frag, int n_read, const char* id, const char* seq, const char* qual)
{
    
}
