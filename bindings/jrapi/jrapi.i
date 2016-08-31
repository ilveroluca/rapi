

/* Java SWIG interface for rapi.h */

%module Rapi

/*
Remove rapi_ prefix from all function and structure names.
Make sure this general rename comes before the ignore clauses below
or they won't have effect.
*/
%rename("NO_ERROR")               "RAPI_NO_ERROR";
%rename("GENERIC_ERROR")          "RAPI_GENERIC_ERROR";
%rename("OP_NOT_SUPPORTED_ERROR") "RAPI_OP_NOT_SUPPORTED_ERROR";
%rename("MEMORY_ERROR")           "RAPI_MEMORY_ERROR";
%rename("PARAM_ERROR")            "RAPI_PARAM_ERROR";
%rename("TYPE_ERROR")             "RAPI_TYPE_ERROR";


%rename("%(strip:[rapi_])s") ""; // e.g., rapi_load_ref -> load_ref


// ignore functions from klib
%rename("$ignore", regextarget=1) "^k";
%rename("$ignore", regextarget=1) "^K";

// also ignore anything that starts with an underscore
%rename("$ignore", regextarget=1) "^_";

// Swig doesn't support applying multiple renaming rules to the same identifier.
// To get camel case, we'll need to rename the wrapped classes individually.
// %rename("%(camelcase)s", %$isclass) "";

// First we need to include the helpers, since they're used throughput
// the wrapping code
%include "../rapi_swig_helpers.i"

/**************************************************
 ###  Typemaps
 **************************************************/

/*
rapi_error_t is the error type for the RAPI library.
*/
typedef int rapi_error_t;
typedef long long rapi_ssize_t;

%inline %{
typedef int rapi_bool;
%}

%typemap(jtype) rapi_bool "boolean";
%typemap(jstype) rapi_bool "boolean";
%typemap(jni) rapi_bool "jboolean";


%mutable;

/**************************************************
 ###  Include rapi_common.i
 * Include it after we've defined all renaming rules, features,
 * typemaps, etc.
 **************************************************/
%include "../rapi_common.i"
/**************************************************
 **************************************************/


%immutable;

/************ begin wrapping functions and structures **************/

// LP: can we avoid redefining the value of the constant here?
#define RAPI_NO_ERROR                    0
#define RAPI_GENERIC_ERROR              -1
#define RAPI_OP_NOT_SUPPORTED_ERROR     -20
#define RAPI_MEMORY_ERROR               -30
#define RAPI_PARAM_ERROR                -40
#define RAPI_TYPE_ERROR                 -50



%mutable;
/********* rapi_opts *******/
typedef struct rapi_opts {
  int ignore_unsupported;
  /* Standard Ones - Differently implemented by aligners*/
  int mapq_min;
  int isize_min;
  int isize_max;
  int n_threads;
  rapi_bool share_ref_mem;

  /* Mismatch / Gap_Opens / Quality Trims --> Generalize ? */

  // TODO: how to wrap this thing?
  //kvec_t(rapi_param) parameters;
} rapi_opts;


///* Init and tear down library */
rapi_error_t rapi_init(const rapi_opts* opts);
rapi_error_t rapi_shutdown(void);


/****** IMMUTABLE *******/
// Everything from here down is read-only

%immutable;

/*
The char* returned by the following functions are wrapped automatically by
SWIG -- the wrapper doesn't try to free the strings since they are "const".
*/
const char* rapi_aligner_name(void);
const char* rapi_aligner_version(void);
const char* rapi_plugin_version(void);

/***************************************
 ****** rapi_contig and rapi_ref *******
 ***************************************/


typedef struct rapi_contig {
  char * name;
  uint32_t len;
  char * assembly_identifier;
  char * species;
  char * uri;
  char * md5;
} rapi_contig;


typedef struct rapi_ref {
  char * path;
  int n_contigs;
} rapi_ref;


/* Load reference */
rapi_error_t rapi_ref_load( const char * reference_path, rapi_ref * ref_struct );

/* Free reference */
rapi_error_t rapi_ref_free( rapi_ref * ref_struct );


/***************************************/
/*      Alignments                     */
/***************************************/

typedef struct rapi_alignment {
  rapi_contig* contig;
  unsigned long int pos; // 1-based
  uint8_t mapq;
  int score; // aligner-specific score

  uint8_t n_mismatches;
  uint8_t n_gap_opens;
  uint8_t n_gap_extensions;
} rapi_alignment;


/***************************************/
/*      Reads and read batches         */
/***************************************/
typedef struct rapi_read {
  char * id;
  char * seq;
  char * qual;
  unsigned int length;
} rapi_read;


typedef struct rapi_batch {
  rapi_ssize_t n_frags;
  int n_reads_frag;
} rapi_batch;

/* Allocate and load reads */
rapi_error_t rapi_reads_alloc( rapi_batch * batch, int n_reads_fragment, int n_fragments );
rapi_error_t rapi_reads_reserve(rapi_batch* batch, rapi_ssize_t n_fragments);
rapi_error_t rapi_reads_clear(rapi_batch* batch);
rapi_error_t rapi_reads_free(rapi_batch * batch);
rapi_ssize_t rapi_batch_read_capacity(const rapi_batch* batch);
rapi_error_t rapi_set_read(rapi_batch * batch, rapi_ssize_t n_frag, int n_read, const char* id, const char* seq, const char* qual, int q_offset);
rapi_read* rapi_get_read(const rapi_batch* batch, rapi_ssize_t n_frag, int n_read);

//rapi_error_t rapi_aligner_state_init(struct rapi_aligner_state** ret_state, const rapi_opts* opts);
//
//rapi_error_t rapi_align_reads( const rapi_ref* ref, rapi_batch* batch,
//    rapi_ssize_t start_frag, rapi_ssize_t end_frag, rapi_aligner_state* state );
//
//rapi_error_t rapi_aligner_state_free(struct rapi_aligner_state* state);
//
///******* SAM output *******/
//rapi_error_t rapi_format_sam(const rapi_read** reads, int n_reads, kstring_t* output);
//rapi_error_t rapi_format_sam_b(const rapi_batch* batch, rapi_ssize_t n_frag, kstring_t* output);
//rapi_error_t rapi_format_sam_hdr(const rapi_ref* ref, kstring_t* output);


/**************************************************
 ###  Object extensions
 **************************************************/

/***************************************
 ****** other stuff              *******
 ***************************************/
