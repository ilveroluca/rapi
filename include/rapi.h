/*
 * RAPI - the Read aligner API
 *
 * Authors:  Luca Pireddu, Riccardo Berutti
 */

#ifndef __RAPI_H__
#define __RAPI_H__

#define RAPI_API_VERSION "0.0"

#include <stdint.h>
#include <kstring.h>
#include <kvec.h>

typedef int rapi_error_t;

/* Error types */
#define RAPI_NO_ERROR                    0
#define RAPI_GENERIC_ERROR              -1
#define RAPI_OP_NOT_SUPPORTED_ERROR     -20
#define RAPI_MEMORY_ERROR               -30
#define RAPI_PARAM_ERROR                -40
#define RAPI_TYPE_ERROR                 -50

static const char*const rapi_error_name(rapi_error_t e)
{
	switch (e) {
		case RAPI_NO_ERROR:               return "NO_ERROR";
		case RAPI_GENERIC_ERROR:          return "GENERIC_ERROR";
		case RAPI_OP_NOT_SUPPORTED_ERROR: return "OP_NOT_SUPPORTED_ERROR";
		case RAPI_MEMORY_ERROR:           return "MEMORY_ERROR";
		case RAPI_PARAM_ERROR:            return "PARAM_ERROR";
		case RAPI_TYPE_ERROR:             return "TYPE_ERROR";
		default: return "Unknown error type";
	};
}

static const char rapi_cigops_char[] = "MIDSH";

/* Key-value TYPES */

#define RAPI_VTYPE_CHAR       1
#define RAPI_VTYPE_TEXT       2
#define RAPI_VTYPE_INT        3
#define RAPI_VTYPE_REAL       4

/* Constants */

#define RAPI_QUALITY_ENCODING_SANGER   33
#define RAPI_QUALITY_ENCODING_ILLUMINA 64
#define RAPI_MAX_TAG_LEN                6

/************************* parameter and tag structures and functions **************/

static inline void rapi_kstr_init(kstring_t* s) {
	s->l = s->m = 0;
	s->s = NULL;
}


typedef struct {
	kstring_t name;
	uint8_t type;
	union {
		char character;
		char* text; // if set, type should be set to TEXT and the str will be freed by rapi_param_clear
		long integer;
		double real;
	} value;
} rapi_param;

static inline void rapi_param_init(rapi_param* kv) {
	memset(kv, 0, sizeof(*kv));
}

static inline void rapi_param_free(rapi_param* kv) {
	free(kv->name.s);
	kv->name.l = kv->name.m = 0;
	kv->name.s = NULL;
	if (kv->type == RAPI_VTYPE_TEXT) {
		free(kv->value.text);
		kv->value.text = NULL;
	}
}

static inline void rapi_param_set_name( rapi_param* kv, const char* key) {
  kv->name.l = 0; // reset cursor position to start
  kputs(key, &kv->name);
}

#define KV_SET_IMPL(value_type, value_field) \
{\
	kv->type = (value_type);\
	(kv->value_field) = value;\
}


static inline void rapi_param_set_char(rapi_param* kv, char value   ) KV_SET_IMPL(RAPI_VTYPE_CHAR, value.character)
// set_text does not take ownership of the string
static inline void rapi_param_set_text(rapi_param* kv, char* value  ) KV_SET_IMPL(RAPI_VTYPE_TEXT, value.text)
static inline void rapi_param_set_long(rapi_param* kv, long value   ) KV_SET_IMPL(RAPI_VTYPE_INT,  value.integer)
static inline void rapi_param_set_dbl( rapi_param* kv, double value ) KV_SET_IMPL(RAPI_VTYPE_REAL, value.real)

static inline const char* rapi_param_get_name(const rapi_param* kv) { return kv->name.s; }

#define KV_GET_IMPL(value_type, value_field) \
{\
	if (kv->type == (value_type)) {\
		*value = (kv->value_field);\
		return RAPI_NO_ERROR;\
	}\
	else\
		return RAPI_TYPE_ERROR;\
}

static inline int rapi_param_get_char(const rapi_param* kv, char * value      ) KV_GET_IMPL(RAPI_VTYPE_CHAR, value.character)
static inline int rapi_param_get_text(const rapi_param* kv, const char** value) KV_GET_IMPL(RAPI_VTYPE_TEXT, value.text)
static inline int rapi_param_get_long(const rapi_param* kv, long * value      ) KV_GET_IMPL(RAPI_VTYPE_INT,  value.integer)
static inline int rapi_param_get_dbl( const rapi_param* kv, double * value    ) KV_GET_IMPL(RAPI_VTYPE_REAL, value.real)

/* Key-value list */
typedef struct {
	char key[RAPI_MAX_TAG_LEN + 1]; // null-terminated
	uint8_t type;
	union {
		char character;
		kstring_t text;
		long integer;
		double real;
	} value;
} rapi_tag;


static inline void rapi_tag_set_key(rapi_tag* kv, const char* s) {
	strncpy(kv->key, s, RAPI_MAX_TAG_LEN);
	kv->key[RAPI_MAX_TAG_LEN] = '\0'; // null terminate, always
}

static inline void rapi_tag_clear(rapi_tag* kv) {
	if (kv->type == RAPI_VTYPE_TEXT) {
		free(kv->value.text.s);
		kv->value.text.s = NULL;
	}
	kv->type = 0;
}

/**
 * Set the tag value to TEXT type and copy `value` into it.
 *
 * If you need to write non-string data, call this fn with a `value` of ""
 * and then use the kstring functions directly with kv->value.text.
 *
 * NOTE: if you set_text and then call any other set_ function without calling
 * rapi_tag_clear you'll leak the previous string value (it won't be automatically
 * cleared).
 */
static inline void rapi_tag_set_text(rapi_tag* kv, const char* value) {
	kv->type = RAPI_VTYPE_TEXT;
	rapi_kstr_init(&kv->value.text);
	kputs(value, &kv->value.text);
}

static inline void rapi_tag_set_char(rapi_tag* kv, char value       ) KV_SET_IMPL(RAPI_VTYPE_CHAR, value.character)
static inline void rapi_tag_set_long(rapi_tag* kv, long value       ) KV_SET_IMPL(RAPI_VTYPE_INT,  value.integer)
static inline void rapi_tag_set_dbl( rapi_tag* kv, double value     ) KV_SET_IMPL(RAPI_VTYPE_REAL, value.real)

static inline int rapi_tag_get_text(const rapi_tag* kv, const kstring_t** value) {
	if (kv->type == RAPI_VTYPE_TEXT) {
		*value = &kv->value.text;
		return RAPI_NO_ERROR;
	}
	else
		return RAPI_TYPE_ERROR;
}

static inline int rapi_tag_get_char(const rapi_tag* kv, char * value      ) KV_GET_IMPL(RAPI_VTYPE_CHAR, value.character)
static inline int rapi_tag_get_long(const rapi_tag* kv, long * value      ) KV_GET_IMPL(RAPI_VTYPE_INT,  value.integer)
static inline int rapi_tag_get_dbl( const rapi_tag* kv, double * value    ) KV_GET_IMPL(RAPI_VTYPE_REAL, value.real)


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
	kvec_t(rapi_param) parameters;

	void * _private; /**< can be used for aligner-specific data */
} rapi_opts;


/************************* The meat starts **************/

/**
 * Reference
 */
typedef struct {
	char * name;
	uint32_t len;
	char * assembly_identifier;
	char * species;
	char * uri;
	char * md5;
} rapi_contig;

typedef struct {
	char * path;
	int n_contigs;
	rapi_contig * contigs;
	void * _private;
} rapi_ref;


/**
 * Alignment
 */
typedef struct {
	uint32_t op:4,
	         len:28;
} rapi_cigar;

typedef kvec_t(rapi_tag) rapi_tag_list;

typedef struct {
	rapi_contig* contig;
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

	rapi_cigar * cigar_ops;
	uint8_t n_cigar_ops;

	rapi_tag_list tags;
} rapi_alignment;

/**
 * Reads
 */
typedef struct {
	char * id;   // NULL-terminated
	char * seq;  // NULL-terminated, capital letters in [AGCTN]
	char * qual; // NULL-terminated, ASCII-encoded in Sanger q+33 format
	unsigned int length; // sequence length
	rapi_alignment* alignments;
	uint8_t n_alignments;
} rapi_read;

/**
 * Batches of reads
 */
typedef struct {
	int n_frags;
	int n_reads_frag;
	rapi_read * reads;
} rapi_batch;


/*************************** functions *******************************/

/* Init Options */
rapi_error_t rapi_opts_init( rapi_opts * my_opts );

rapi_error_t rapi_opts_free( rapi_opts * my_opts );

/* Init and tear down library */
rapi_error_t rapi_init(const rapi_opts* opts);
rapi_error_t rapi_shutdown();

/* Aligner Version */
const char* rapi_aligner_name();
const char* rapi_aligner_version();
const char* rapi_plugin_version();

/* Load reference */
rapi_error_t rapi_ref_load( const char * reference_path, rapi_ref * ref_struct );

/* Free reference */
rapi_error_t rapi_ref_free( rapi_ref * ref_struct );

/* Allocate reads */
rapi_error_t rapi_reads_alloc( rapi_batch * batch, int n_reads_fragment, int n_fragments );

/* Reserve sufficient space for n_fragments.
 *
 * If n_fragments is greater than the number of fragments that can be currently
 * stored in the rapi_batch (i.e., greater than n_reads_frag*n_frags), this
 * function reallocates the space to fit them.  The new space will be
 * initialized to 0, like rapi_reads_alloc * does.
 *
 * This function will not reduce the space occupied by the batch;
 * if n_fragments is less than what can be already stored then the function
 * does nothing.
 *
 * In case of error, the batch is not modified.
 */
rapi_error_t rapi_reads_reserve(rapi_batch* batch, int n_fragments);

/**
 * Empty a `batch`, clearing any reads stored therein.  The `batch` will be
 * restored to a state as if it was just allocated by rapi_reads_alloc
 * (so the read memory is not freed).
 */
rapi_error_t rapi_reads_clear(rapi_batch* batch);

/**
 * Clear a read `batch` and free all associated memory.
 */
rapi_error_t rapi_reads_free( rapi_batch * batch );


/**
 * Number of reads that fit in curretly allocated space.
 */
static inline int rapi_batch_read_capacity(const rapi_batch* batch) {
	return batch->n_frags * batch->n_reads_frag;
}

/**
 * Set read data within a batch.  The strings are copied into the read batch.
 *
 * \param n_frag 0-based fragment number
 * \param n_read 0-based read number
 * \param id read name (NULL-terminated)
 * \param seq base sequence (NULL-terminated)
 * \param qual per-base quality, or NULL
 * \param q_offset offset from 0 for the base quality values (e.g., 33 for Sanger, 0 for byte values)
 */
rapi_error_t rapi_set_read(rapi_batch * batch, int n_frag, int n_read, const char* id, const char* seq, const char* qual, int q_offset);


/* Align */
typedef struct rapi_aligner_state rapi_aligner_state; //< opaque structure.  Aligner can use for whatever it wants.

rapi_error_t rapi_aligner_state_init(const rapi_opts* opts, struct rapi_aligner_state** ret_state);

/**
 * Align the reads in batch to ref.
 *
 * The alignments are written directly to the read structures in the rapi_batch.
 *
 * \param ref The alignment reference
 * \param batch A read batch
 * \param start_frag Within the read batch, start aligning reads at this fragment (0-based)
 * \param end_frag Stop aligning at this fragment (exclusive). So, to align only the
 *                 second pair of reads in the batch give the indices [1, 2). For the entire
 *                 batch give [0, batch.n_frags).
 * \param state Provide the state initialized with rapi_aligner_state_init.
 */
rapi_error_t rapi_align_reads( const rapi_ref* ref, rapi_batch* batch,
    int start_frag, int end_frag, rapi_aligner_state* state );

rapi_error_t rapi_aligner_state_free(struct rapi_aligner_state* state);

static inline rapi_read* rapi_get_read(const rapi_batch* batch, int n_frag, int n_read) {
	if (n_frag >= 0 || n_frag < batch->n_frags
	 || n_read >= 0 || n_read < batch->n_reads_frag) {
	  return batch->reads + (n_frag * batch->n_reads_frag + n_read);
  }
  return NULL; // else coordinates are out of bounds
}

long rapi_get_insert_size(const rapi_alignment* read, const rapi_alignment* mate);

int rapi_get_rlen(int n_cigar, const rapi_cigar* cigar_ops);

void rapi_put_cigar(int n_ops, const rapi_cigar* ops, int force_hard_clip, kstring_t* output);


/******* SAM output *******/

/**
 * Format SAM for all reads in the given fragment.
 *
 * \param batch A read batch.
 *
 * \param n_frag The index of the fragment within the batch whose reads will be
 *               printed in the SAM text.  SAM will be produced for all the reads
 *               within the fragment.
 *
 * \param output An initialized kstring_t to which the SAM will be appended.
 */
rapi_error_t rapi_format_sam(const rapi_batch* batch, int n_frag, kstring_t* output);

/**
 * Format the SAM header for the given reference.  The header will also contain
 * a @PG tag identifying the RAPI-interfaced aligner being used.
 *
 * \param output An initialized kstring_t to which the output will be appended.
 */
rapi_error_t rapi_format_sam_hdr(const rapi_ref* ref, kstring_t* output);

#endif
