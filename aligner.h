/*
  aligner.h
*/

#define ALN_TAG_TYPE_CHAR		1
#define ALN_TAG_TYPE_TEXT		2
#define ALN_TAG_TYPE_INT		3
#define ALN_TAG_TYPE_REAL		4
#define ALN_TAG_TYPE_BYTE_ARRAY		5
#define ALN_TAG_TYPE_INT_ARRAY		6

#define ALN_QUALITY_ENCODING_SANGER	33
#define ALN_QUALITY_ENCODING_ILLUMINA	64

/* Structs */
  typedef struct {
    int ignore_unsupported;
    /* Standard Ones - Differently implemented by aligners*/
    uint8_t mapq_min;
    uint8_t trim_quality;
    uint8_t isize_min;
    uint16_t isize_max;
    uint8_t seed_len;
    /* Mismatch / Gap_Opens / Quality Trims --> Generalize ? */
    /* Aligner specific parameters */
    int n_parameters;
    aln_parameter * parameters;
  } aln_config;
  
  typedef struct {
     char * name;
     void * value;
  } aln_parameter
  } 

  typedef struct {
    const char * path;
    int n_contigs;
    aln_contig * contigs;
    void * my_reference;
  } aln_ref;

  typedef struct {
    const char * name;
    uint32_t length;
    const char * assembly_identifier;
    const char md5[32];
    const char * species;
    const char * uri;
  } aln_contig;
    
  typedef struct {
    int n_allocated_frags;
    int n_frags;
    int n_reads_frag;
    aln_read * reads;
  } aln_batch;


  typedef struct {
    char * id;
    unsigned char * seq;
    unsigned char * qual;
    unsigned int length;
    unsigned int bp_clipping;
    aln_alignment alignment;
  } aln_read;


  typedef struct {
    char * contig;
    unsigned long int pos;
    aln_tag * tags;
    uint8_t n_tags;
    uint8_t mapq;
    aln_cigar * cigar_ops;
    uint8_t n_cigar_ops;
    uint8_t	paired:1,
		prop_paired:1,
		unmapped:1,
		reverse_strand:1,  
		secondary_aln:1;
    uint8_t	n_mismatches;
    uint8_t	n_gap_opens;
    uint8_t	n_gap_extensions;
    uint16_t	edit_distance;
    aln_mismatch * mismatches;
    uint8_t	n_mismatches;
  } aln_alignment;

  
  typedef struct {
    uint32_t op:4,
	     len:28;
  } aln_cigar;
  
  typedef struct {
    uint8_t type;
    uint8_t len;
    char * data;
  } aln mismatch;
  
  typedef struct {
    char * name;
    unsigned byte type;
    union {
		  char character,
		  char * text,
		  long integer,
		  double real,
		  uint8_t * byte_array,
		  long * integer_array
    } data;
  } aln_tag


/* Init Library */
int aln_init();

/* Init Options */
int aln_init_opts( aln_opts * my_opts );

/* Aligner Version */
const char * aln_version();

/* Load reference */
int aln_load_ref( const char * reference, aln_ref * ref_struct );

/* Free reference */
int aln_free_ref( aln_ref * ref_struct );

/* Allocate reads */
int aln_alloc_reads( aln_batch ** batch, int n_reads_fragment, int n_fragments );

/* Free reads */
int aln_free_reads( aln_batch ** batch );

/* Align */
int aln_align_reads( aln_batch * batch, aln_config * config );

