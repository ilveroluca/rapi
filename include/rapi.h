/*
 * RAPI - the Read aligner API
 */

/******************************************************************************
 *  Copyright (c) 2014-2016 Center for Advanced Studies,
 *                          Research and Development in Sardinia (CRS4)
 *  
 *  Licensed under the terms of the MIT License (see LICENSE file included with the
 *  project).
 *  
 *  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 *  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 *  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 *  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 *  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 *  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 *  SOFTWARE.
 *****************************************************************************/

#ifndef __RAPI_H__
#define __RAPI_H__

#define RAPI_API_VERSION "0.0"

#include <stddef.h>
#include <kstring.h>
#include <kvec.h>

typedef long long rapi_ssize_t;

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

#define RAPI_CIG_M            0 // Match
#define RAPI_CIG_I            1 // Insertion
#define RAPI_CIG_D            2 // Deletion
#define RAPI_CIG_S            3 // Soft-clip
#define RAPI_CIG_H            4 // Hard-clip
#define RAPI_CIG_N            5 // Skipped region from the reference
#define RAPI_CIG_P            6 // Padding

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


typedef struct rapi_param {
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
typedef struct rapi_tag {
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
typedef struct rapi_opts {

  /** Tell implementation to ignore unsupported options.
   * Alternatively, it should give an error */
	int ignore_unsupported;
	// alignment filtering
	int mapq_min;
	int isize_min;
	int isize_max;

	// multithreading -- implementation may ignore it if single-threaded
	int n_threads;

	// Whether to share references in memory with other processes using RAPI
	// (if the implementation supports it)
	int share_ref_mem;

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
typedef struct rapi_contig {
	char * name;
	rapi_ssize_t len;
	char * assembly_identifier;
	char * species;
	char * uri;
	char * md5;
} rapi_contig;

typedef struct rapi_ref {
	char * path;
	int n_contigs;
	rapi_contig * contigs;
	void * _private;
} rapi_ref;


/**
 * Alignment
 */
typedef struct rapi_cigar {
	uint32_t op:4,
	         len:28;
} rapi_cigar;

typedef kvec_t(rapi_tag) rapi_tag_list;

typedef struct rapi_alignment {
	rapi_contig* contig;
	rapi_ssize_t pos; // 1-based
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
typedef struct rapi_read {
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
typedef struct rapi_batch {
	rapi_ssize_t n_frags;
	int n_reads_frag;
	void * _private;
} rapi_batch;


/*************************** functions *******************************/

/** Initialize rapi_opts structure */
rapi_error_t rapi_opts_init( rapi_opts * my_opts );
/** Destroy rapi_opts structure */
rapi_error_t rapi_opts_free( rapi_opts * my_opts );

/** Initialize library.
 *
 * \param opts Options may be provided to the implementation.
 *
 * If options are not provided (i.e, `opts` is NULL), the library will use
 * default values.
 *
 * If options are provided, they may be remembered by the library as needed in
 * other calls.
 */
rapi_error_t rapi_init(const rapi_opts* opts);

/** Shutdown library, clearing any allocated resources. */
rapi_error_t rapi_shutdown(void);

/** Get wrapped aligner name */
const char* rapi_aligner_name(void);
/** Get wrapped aligner version */
const char* rapi_aligner_version(void);
/** Get aligner plug-in version */
const char* rapi_plugin_version(void);

/** Load a reference.
 *
 * The implementation may configure its behaviour based on the options passed
 * into rapi_init.
 */
rapi_error_t rapi_ref_load( const char * reference_path, rapi_ref * ref_struct );

/** Free reference structure and unload reference (if loaded). */
rapi_error_t rapi_ref_free( rapi_ref * ref_struct );

/**
 * Create read batch configured for `n_reads_fragment` reads per fragment.
 * Allocate memory for `n_fragments` fragments.
 */
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
rapi_error_t rapi_reads_reserve(rapi_batch* batch, rapi_ssize_t n_fragments);

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
static inline rapi_ssize_t rapi_batch_read_capacity(const rapi_batch* batch) {
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
rapi_error_t rapi_set_read(rapi_batch * batch, rapi_ssize_t n_frag, int n_read, const char* id, const char* seq, const char* qual, int q_offset);

/**
 * Get pointer to read at coordinates (n_frag, n_read).
 *
 * \param n_frag Index in interval [0, batch->n_frags).
 * \param n_read Index in interval [0, batch->n_reads_frag).
 *
 * \return pointer to rapi_read structure;  NULL in case of error.
 *
 * \warning The returned pointer may point to an uninitialized structure if the
 * space has been allocated (with rapi_reads_reserve or rapi_reads_alloc) but
 * not set with rapi_set_read.
 */
rapi_read* rapi_get_read(const rapi_batch* batch, rapi_ssize_t n_frag, int n_read);


/* Aligner section */

/** Opaque aligner structure.  Aligner can use for whatever it wants. */
typedef struct rapi_aligner_state rapi_aligner_state;

/**
 * Create the aligner state structure.
 *
 * The library should initialize the aligner based on the options passed to
 * rapi_init().  However, the user can provide new options that override the
 * library-wide configuration.  Alternatively, \param opts is NULL.
 *
 * \param ret_state Return argument for the new aligner state.
 * \param opts User-specified options, or NULL if the user wants to use the options passed to rapi_init.
 * \return rapi_error_t Return code.
 */
rapi_error_t rapi_aligner_state_init(struct rapi_aligner_state** ret_state, const rapi_opts* opts);

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
    rapi_ssize_t start_frag, rapi_ssize_t end_frag, rapi_aligner_state* state );

/** Clear aligner state and free any associated system resources. */
rapi_error_t rapi_aligner_state_free(struct rapi_aligner_state* state);

#endif
