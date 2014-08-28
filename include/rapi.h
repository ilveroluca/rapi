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

/* Error types */
#define RAPI_NO_ERROR                    0
#define RAPI_GENERIC_ERROR              -1
#define RAPI_OP_NOT_SUPPORTED_ERROR     -20
#define RAPI_MEMORY_ERROR               -30
#define RAPI_PARAM_ERROR                -40
#define RAPI_TYPE_ERROR                 -50

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


static inline void rapi_param_set_char(rapi_param* kv, char value       ) KV_SET_IMPL(RAPI_VTYPE_CHAR, value.character)
static inline void rapi_param_set_text(rapi_param* kv, char* value) KV_SET_IMPL(RAPI_VTYPE_TEXT, value.text)
static inline void rapi_param_set_long(rapi_param* kv, long value       ) KV_SET_IMPL(RAPI_VTYPE_INT,  value.integer)
static inline void rapi_param_set_dbl( rapi_param* kv, double value     ) KV_SET_IMPL(RAPI_VTYPE_REAL, value.real)

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

	kvec_t(rapi_tag) tags;
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
int rapi_opts_init( rapi_opts * my_opts );

int rapi_opts_free( rapi_opts * my_opts );

/* Init and tear down library */
int rapi_init(const rapi_opts* opts);
int rapi_shutdown();

/* Aligner Version */
const char* rapi_aligner_name();
const char* rapi_aligner_version();

/* Load reference */
int rapi_ref_load( const char * reference_path, rapi_ref * ref_struct );

/* Free reference */
int rapi_ref_free( rapi_ref * ref_struct );

/* Allocate reads */
int rapi_reads_alloc( rapi_batch * batch, int n_reads_fragment, int n_fragments );

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
int rapi_reads_reserve(rapi_batch * batch, int n_fragments);

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
int rapi_set_read(rapi_batch * batch, int n_frag, int n_read, const char* id, const char* seq, const char* qual, int q_offset);

int rapi_reads_free( rapi_batch * batch );

/* Align */
typedef struct rapi_aligner_state rapi_aligner_state; //< opaque structure.  Aligner can use for whatever it wants.

int rapi_aligner_state_init(const rapi_opts* opts, struct rapi_aligner_state** ret_state);

int rapi_align_reads( const rapi_ref* ref,  rapi_batch * batch, const rapi_opts * config, rapi_aligner_state* state );

int rapi_aligner_state_free(struct rapi_aligner_state* state);

static inline rapi_read* rapi_get_read(const rapi_batch* batch, int n_fragment, int n_read) {
	return batch->reads + (n_fragment * batch->n_reads_frag + n_read);
}

long rapi_get_insert_size(const rapi_alignment* read, const rapi_alignment* mate);

static inline int rapi_get_rlen(int n_cigar, const rapi_cigar* cigar_ops)
{
	int len = 0;
	for (int k = 0; k < n_cigar; ++k) {
		int op = cigar_ops[k].op;
		if (op == 0 || op == 2)
			len += cigar_ops[k].len;
	}
	return len;
}

int rapi_format_sam(const rapi_read* read, const rapi_read* mate, kstring_t* output);

void rapi_put_cigar(int n_ops, const rapi_cigar* ops, int force_hard_clip, kstring_t* output);

#endif
