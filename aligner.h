/*
 * RAPI - the Read aligner API
 *
 * Authors:  Luca Pireddu, Riccardo Berutti, Ridvan Dongelci
 */

#ifndef __RAPI_H__
#define __RAPI_H__

#include <stdint.h>
#include <kstring.h>
#include <kvec.h>

/* Error types */
#define ALN_NO_ERROR                     0
#define ALN_GENERIC_ERROR               -1
#define ALN_OP_NOT_SUPPORTED_ERROR      -2

#define ALN_REFERENCE_ERROR            -10

#define ALN_TAG_NOT_EXISTING           -20

#define ALN_MEMORY_ERROR               -30
#define ALN_PARAM_ERROR                -40
#define ALN_TYPE_ERROR                 -50

/* Key-value TYPES */

#define ALN_VTYPE_CHAR       1
#define ALN_VTYPE_TEXT       2
#define ALN_VTYPE_INT        3
#define ALN_VTYPE_REAL       4

/* Constants */

#define ALN_QUALITY_ENCODING_SANGER   33
#define ALN_QUALITY_ENCODING_ILLUMINA 64
#define ALN_MAX_TAG_LEN                6


static inline void aln_init_kstr(kstring_t* s) {
	s->l = s->m = 0;
	s->s = NULL;
}


typedef struct aln_param {
	kstring_t name;
	uint8_t type;
	union {
		char character;
		char* text; // if set, type should be set to TEXT and the str will be freed by aln_param_clear
		long integer;
		double real;
	} value;
} aln_param;

static inline void aln_param_init(aln_param* kv) {
	memset(kv, 0, sizeof(*kv));
}

static inline void aln_param_clear(aln_param* kv) {
	free(kv->name.s);
	kv->name.l = kv->name.m = 0;
	kv->name.s = NULL;
	if (kv->type == ALN_VTYPE_TEXT) {
		free(kv->value.text);
		kv->value.text = NULL;
	}
}

static inline void aln_param_set_name( aln_param* kv, const char* key) { kputs(key, &kv->name); }

#define KV_SET_IMPL(value_type, value_field) \
{\
	kv->type = (value_type);\
	(kv->value_field) = value;\
}


static inline void aln_param_set_char(aln_param* kv, char value       ) KV_SET_IMPL(ALN_VTYPE_CHAR, value.character)
static inline void aln_param_set_text(aln_param* kv, char* value) KV_SET_IMPL(ALN_VTYPE_TEXT, value.text)
static inline void aln_param_set_long(aln_param* kv, long value       ) KV_SET_IMPL(ALN_VTYPE_INT,  value.integer)
static inline void aln_param_set_dbl( aln_param* kv, double value     ) KV_SET_IMPL(ALN_VTYPE_REAL, value.real)

static inline const char* aln_param_get_name(const aln_param* kv) { return kv->name.s; }

#define KV_GET_IMPL(value_type, value_field) \
{\
	if (kv->type == (value_type)) {\
		*value = (kv->value_field);\
		return ALN_NO_ERROR;\
	}\
	else\
		return ALN_TYPE_ERROR;\
}

static inline int aln_param_get_char(const aln_param* kv, char * value      ) KV_GET_IMPL(ALN_VTYPE_CHAR, value.character)
static inline int aln_param_get_text(const aln_param* kv, const char** value) KV_GET_IMPL(ALN_VTYPE_TEXT, value.text)
static inline int aln_param_get_long(const aln_param* kv, long * value      ) KV_GET_IMPL(ALN_VTYPE_INT,  value.integer)
static inline int aln_param_get_dbl( const aln_param* kv, double * value    ) KV_GET_IMPL(ALN_VTYPE_REAL, value.real)

/* Key-value list */
typedef struct aln_tag {
	char key[ALN_MAX_TAG_LEN + 1]; // null-terminated
	uint8_t type;
	union {
		char character;
		kstring_t text;
		long integer;
		double real;
	} value;
} aln_tag;


static inline void aln_tag_set_key(aln_tag* kv, const char* s) {
	strncpy(kv->key, s, ALN_MAX_TAG_LEN);
	kv->key[ALN_MAX_TAG_LEN] = '\0'; // null terminate, always
}

static inline void aln_tag_clear(aln_tag* kv) {
	if (kv->type == ALN_VTYPE_TEXT) {
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
static inline void aln_tag_set_text(aln_tag* kv, const char* value) {
	kv->type = ALN_VTYPE_TEXT;
	aln_init_kstr(&kv->value.text);
	kputs(value, &kv->value.text);
}

static inline void aln_tag_set_char(aln_tag* kv, char value       ) KV_SET_IMPL(ALN_VTYPE_CHAR, value.character)
static inline void aln_tag_set_long(aln_tag* kv, long value       ) KV_SET_IMPL(ALN_VTYPE_INT,  value.integer)
static inline void aln_tag_set_dbl( aln_tag* kv, double value     ) KV_SET_IMPL(ALN_VTYPE_REAL, value.real)

static inline int aln_tag_get_text(const aln_tag* kv, const kstring_t** value) {
	if (kv->type == ALN_VTYPE_TEXT) {
		*value = &kv->value.text;
		return ALN_NO_ERROR;
	}
	else
		return ALN_TYPE_ERROR;
}

static inline int aln_tag_get_char(const aln_tag* kv, char * value      ) KV_GET_IMPL(ALN_VTYPE_CHAR, value.character)
static inline int aln_tag_get_long(const aln_tag* kv, long * value      ) KV_GET_IMPL(ALN_VTYPE_INT,  value.integer)
static inline int aln_tag_get_dbl( const aln_tag* kv, double * value    ) KV_GET_IMPL(ALN_VTYPE_REAL, value.real)


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
	kvec_t(aln_param) parameters;

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

	aln_cigar * cigar_ops;
	uint8_t n_cigar_ops;

	kvec_t(aln_tag) tags;
} aln_alignment;

typedef struct {
	char * id;   // NULL-terminated
	char * seq;  // NULL-terminated, capital letters in [AGCTN]
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

int aln_free_reads( aln_batch * batch );

/* Align */
typedef struct aln_aligner_state aln_aligner_state; //< opaque structure.  Aligner can use for whatever it wants.

int aln_init_aligner_state(const aln_opts* opts, struct aln_aligner_state** ret_state);

int aln_align_reads( const aln_ref* ref,  aln_batch * batch, const aln_opts * config, aln_aligner_state* state );

int aln_free_aligner_state(struct aln_aligner_state* state);

static inline aln_read* aln_get_read(const aln_batch* batch, int n_fragment, int n_read) {
	return batch->reads + (n_fragment * batch->n_reads_frag + n_read);
}

long aln_get_insert_size(const aln_alignment* read, const aln_alignment* mate);

static inline int aln_get_rlen(int n_cigar, const aln_cigar* cigar_ops)
{
	int len = 0;
	for (int k = 0; k < n_cigar; ++k) {
		int op = cigar_ops[k].op;
		if (op == 0 || op == 2)
			len += cigar_ops[k].len;
	}
	return len;
}

int aln_format_sam(const aln_read* read, const aln_read* mate, kstring_t* output);

void aln_put_cigar(int n_ops, const aln_cigar* ops, int force_hard_clip, kstring_t* output);

#endif
