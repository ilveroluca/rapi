/*
 * RAPI - the Read aligner API
 *
 * Authors:  Luca Pireddu, Riccardo Berutti, Ridvan Dongelci
 */

#ifndef __RAPI_H__
#define __RAPI_H__

#include <stdint.h>

/* Error types */
#define ALN_NO_ERROR                     0
#define ALN_GENERIC_ERROR               -1
#define ALN_OP_NOT_SUPPORTED_ERROR      -2

#define ALN_REFERENCE_ERROR            -10

#define ALN_TAG_NOT_EXISTING           -20

#define ALN_MEMORY_ERROR               -30
#define ALN_PARAM_ERROR                -40

/* Key-value TYPES */

#define ALN_VTYPE_CHAR       1
#define ALN_VTYPE_TEXT       2
#define ALN_VTYPE_INT        3
#define ALN_VTYPE_REAL       4
#define ALN_VTYPE_BYTE_ARRAY 5
#define ALN_VTYPE_INT_ARRAY  6

/* FASTQ Quality Encoding */

#define ALN_QUALITY_ENCODING_SANGER 33
#define ALN_QUALITY_ENCODING_ILLUMINA 64


/* Key-value list */
typedef struct aln_kv {
	const char * key;
	uint8_t type;
	union {
		char character;
		const char * text;
		long integer;
		double real;
		const uint8_t * byte_array;
		const long * integer_array;
	} value;

	const struct aln_kv * next;
} aln_kv;

aln_kv* aln_kv_insert_char(const char* key, char value,           const aln_kv* tail);
aln_kv* aln_kv_insert_str( const char* key, const char* value,    const aln_kv* tail);
aln_kv* aln_kv_insert_long(const char* key, long value,           const aln_kv* tail);
aln_kv* aln_kv_insert_dbl( const char* key, double value,         const aln_kv* tail);
aln_kv* aln_kv_insert_ba(  const char* key, const uint8_t* value, const aln_kv* tail);
aln_kv* aln_kv_insert_la(  const char* key, const long* value,    const aln_kv* tail);

int aln_kv_get_char(const aln_kv* kv_list, const char* key, char * value         );
int aln_kv_get_str( const aln_kv* kv_list, const char* key, const char** value   );
int aln_kv_get_long(const aln_kv* kv_list, const char* key, long * value         );
int aln_kv_get_dbl( const aln_kv* kv_list, const char* key, double * value       );
int aln_kv_get_ba(  const aln_kv* kv_list, const char* key, const uint8_t** value);
int aln_kv_get_la(  const aln_kv* kv_list, const char* key, const long** value   );


// free the key and values??
int aln_kv_free_list(aln_kv* list);


/**
 * Options.
 */
typedef struct {
	int ignore_unsupported;
	/* Standard Ones - Differently implemented by aligners*/
	int mapq_min;
	int isize_min;
	int isize_max;
	/* Mismatch / Gap_Opens / Quality Trims --> Generalize ? */

	/* Aligner specific parameters in 'parameters' list.
	 * LP: I'm thinking we might want to drop this list in favour
	 * of letting the user set aligner-specific options through the
	 * _private structure below.
	 */
	int n_parameters;
	aln_kv * parameters;

	void * _private; /**< can be used for aligner-specific data */
} aln_opts;


/**
 * Reference
 */
typedef struct {
	char * name;
	uint32_t len;
	char * assembly_identifier;
	char * species;
	char * uri;
	char md5[32];
} aln_contig;

typedef struct {
	char * path;
	int n_contigs;
	aln_contig * contigs;
	void * _private;
} aln_ref;


/**
 * Read and alignment
 */
typedef struct {
	uint32_t op:4,
	         len:28;
} aln_cigar;

typedef struct {
	uint8_t type;
	uint8_t len;
	char * data;
} aln_mismatch;


typedef struct {
	aln_contig* contig;
	unsigned long int pos; // 1-based
	uint8_t mapq;
	int score; // aligner-specific score

	uint8_t paired:1,
	        prop_paired:1,
	        mapped:1,
	        reverse_strand:1,
	        secondary_aln:1;

	uint8_t n_mismatches;
	uint8_t n_gap_opens;
	uint8_t n_gap_extensions;

	aln_kv * tags;
	aln_cigar * cigar_ops;
	uint8_t n_cigar_ops;
	aln_mismatch * mm_ops;
	uint8_t n_mm_ops;
} aln_alignment;

typedef struct {
	char * id;   // NULL-terminated
	char * seq;  // NULL-terminated
	char * qual; // NULL-terminated, ASCII-encoded in Sanger q+33 format
	unsigned int length;
	aln_alignment* alignments;
	uint8_t n_alignments;
} aln_read;

typedef struct {
	int n_frags;
	int n_reads_frag;
	aln_read * reads;
} aln_batch;

/* Init Options */
int aln_init_opts( aln_opts * my_opts );

int aln_free_opts( aln_opts * my_opts );

/* Init Library */
int aln_init(const aln_opts* opts);

/* Aligner Version */
const char * aln_version();

/* Load reference */
int aln_load_ref( const char * reference_path, aln_ref * ref_struct );

/* Free reference */
int aln_free_ref( aln_ref * ref_struct );

/* Allocate reads */
int aln_alloc_reads( aln_batch * batch, int n_reads_fragment, int n_fragments );

/**
 * Set read data within a batch.
 *
 * \param n_frag 0-based fragment number
 * \param n_read 0-based read number
 * \param id read name (NULL-terminated)
 * \param seq base sequence (NULL-terminated)
 * \param qual per-base quality, or NULL
 * \param q_offset offset from 0 for the base quality values (e.g., 33 for Sanger, 0 for byte values)
 */
int aln_set_read(aln_batch * batch, int n_frag, int n_read, const char* id, const char* seq, const char* qual, int q_offset);

/* Free reads */
int aln_free_reads( aln_batch * batch );

/* Align */
typedef struct aln_aligner_state aln_aligner_state; //< opaque structure.  Aligner can use for whatever it wants.

int aln_init_aligner_state(const aln_opts* opts, struct aln_aligner_state** ret_state);

int aln_align_reads( const aln_ref* ref,  aln_batch * batch, const aln_opts * config, aln_aligner_state* state );

int aln_free_aligner_state(struct aln_aligner_state* state);


inline aln_read* aln_get_read(const aln_batch* batch, int n_fragment, int n_read) {
	return batch->reads + (n_fragment * batch->n_reads_frag + n_read);
}

#endif
