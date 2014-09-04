

/* SWIG interface for rapi.h */

%include "exception.i"

/******* Language-independent exceptions *******/
//SWIG_exception(SWIG_MemoryError, "Not enough memory");
// LP: I can't get SWIG_exception to work. It generates a "goto fail;"
// statement without defining the "fail" label, thus resulting in a
// compiler error.  I've tried this with SWIG 2.0.11 and 3.0.2.
// This leaves me writing one interface per target language.

/******* swig -builtin option *******
PyErr_SetString wouldn't raise a Python exception until I
started calling swig with the "-builtin" option.

*************************************/


%module rapi

%{
/* Includes the header in the wrapper code */
#include "rapi.h"
%}

/*
rapi_error_t is the error type for the RAPI library.
Rather than passing the error codes back to the Python-side caller,
we apply a typemap that:

  * checks whether the code represents an error and, if so, maps it to an
    exception
  * if there's no error, maps the code to None so that the calls don't
    return a value.
*/
typedef int rapi_error_t;

%typemap(out) rapi_error_t {
  if ($1 != RAPI_NO_ERROR) {
    SWIG_exception(rapi_swig_error_type($1), "");
    //PyErr_SetString(rapi_py_error_type($1), "");
  }
  else {
    $result = VOID_Object;
  }
}


%header %{

#include <stddef.h>
#include <stdio.h>

#define PDEBUG(...) { fprintf(stderr, "%s(%d) DEBUG: ", __FILE__, __LINE__); fprintf(stderr, __VA_ARGS__); }
#define PERROR(...) { fprintf(stderr, "%s(%d) ERROR: ", __FILE__, __LINE__); fprintf(stderr, __VA_ARGS__); }

/* Helper to convert RAPI errors to swig errors */
int rapi_swig_error_type(rapi_error_t rapi_code) {
  int type = 0;
  switch(rapi_code) {
  case RAPI_GENERIC_ERROR:
    type = SWIG_RuntimeError;
    break;
  case RAPI_OP_NOT_SUPPORTED_ERROR:
    type = SWIG_RuntimeError;
    break;
  case RAPI_MEMORY_ERROR:
    type = SWIG_MemoryError;
    break;
  case RAPI_PARAM_ERROR:
    type = SWIG_ValueError;
    break;
  case RAPI_TYPE_ERROR:
    type = SWIG_TypeError;
    break;
  default:
    type = SWIG_RuntimeError;
    break;
  };
  return type;
}


/* Helper to convert RAPI errors to Python errors */
PyObject* rapi_py_error_type(rapi_error_t rapi_code) {
  PyObject* type = 0;
  switch(rapi_code) {
  case RAPI_GENERIC_ERROR:
    type = PyExc_RuntimeError;
    break;
  case RAPI_OP_NOT_SUPPORTED_ERROR:
    type = PyExc_NotImplementedError;
    break;
  case RAPI_MEMORY_ERROR:
    type = PyExc_MemoryError;
    break;
  case RAPI_PARAM_ERROR:
    type = PyExc_ValueError;
    break;
  case RAPI_TYPE_ERROR:
    type = PyExc_TypeError;
    break;
  default:
    type = PyExc_RuntimeError;
    break;
  };
  return type;
}

void* rapi_malloc(size_t nbytes) {
  void* result = malloc(nbytes);
  if (!result)
    PyErr_SetString(PyExc_MemoryError, "Failed to allocate memory");

  return result;
}

%}

// We don't provide direct access to the read_batch struct
%rename("$ignore") rapi_batch;
// Instead, we'll define a wrapper for it.  It's implemented
// as the structure rapi_batch_wrap within this wrapper and
// exposed to the user as rapi.read_batch
%rename(read_batch) rapi_batch_wrap;

/*
Remove rapi_ prefix from all function and structure names.
Make sure this general rename comes before the ignore clauses below
or they won't have effect.
*/
%rename("%(strip:[rapi_])s") ""; // e.g., rapi_load_ref -> load_ref


// ignore functions from klib
%rename("$ignore", regextarget=1) "^k";
%rename("$ignore", regextarget=1) "^K";

// also ignore anything that starts with an underscore
%rename("$ignore", regextarget=1) "^_";

%include "stdint.i"; // needed to tell swig about types such as uint32_t, uint8_t, etc.
// This %includes are needed since we have some structure elements that are kvec_t
%include "kvec.h";


// a couple of constants
#define QENC_SANGER   33
#define QENC_ILLUMINA 64

/************ begin functions and structures **************/

rapi_error_t rapi_init(const rapi_opts* opts);
rapi_error_t rapi_shutdown();


/*
The char* returned by the following functions are wrapped automatically by
SWIG -- the wrapper doesn't try to free the strings since they are "const".
*/
const char* rapi_aligner_name();
const char* rapi_aligner_version();

/**** rapi_param ****/
typedef struct {
  uint8_t type;
  union {
    char character;
    char* text; // if set, type should be set to TEXT and the str will be freed by rapi_param_clear
    long integer;
    double real;
  } value;
} rapi_param;


%extend rapi_param {
  const char* name; // kstring_t accessed as a const char*

  rapi_param() {
    rapi_param* param = (rapi_param*) rapi_malloc(sizeof(rapi_param));
    if (!param) return NULL;
    rapi_param_init(param);
    return param;
  }

  ~rapi_param() {
    rapi_param_free($self);
  }
};

%{
const char* rapi_param_name_get(const rapi_param *p) {
  return rapi_param_get_name(p);
}

void rapi_param_name_set(rapi_param *p, const char *val) {
  rapi_param_set_name(p, val);
}
%}

/********* rapi_opts *******/
typedef struct {
  int ignore_unsupported;
  /* Standard Ones - Differently implemented by aligners*/
  int mapq_min;
  int isize_min;
  int isize_max;
  /* Mismatch / Gap_Opens / Quality Trims --> Generalize ? */

  // TODO: how to wrap this thing?
  kvec_t(rapi_param) parameters;
} rapi_opts;


%extend rapi_opts {
  rapi_opts() {
    rapi_opts* opts = (rapi_opts*) rapi_malloc(sizeof(rapi_opts));
    if (!opts) return NULL;

    int error = rapi_opts_init(opts);
    if (error == RAPI_NO_ERROR)
      return opts;
    else {
      free(opts);
      PyErr_SetString(rapi_py_error_type(error), "Not enough memory to initialize rapi_opts");
      return NULL;
    }
  }

  ~rapi_opts() {
    int error = rapi_opts_free($self);
    if (error != RAPI_NO_ERROR) {
      PERROR("Problem destroying opts (error code %d)\n", error);
      // TODO: should we raise exceptions in case of errors when freeing/destroying?
    }
  }
};

/****** rapi_contig and rapi_ref *******/

%immutable; /** don't let the user modify the reference through the wrapper */

typedef struct {
  char * name;
  uint32_t len;
  char * assembly_identifier;
  char * species;
  char * uri;
  char * md5;
} rapi_contig;

/* An iterator object for the contig */
%inline %{
typedef struct {
  rapi_contig* next_item;
  size_t n_left;
} ref_contig_iter;
%}

/* Without these %feature statements the __len__ and __getitem__ methods
   will be defined but they won't work (at least when calling swig with -builtin).
   You'll get errors like:

----> len(r)
TypeError: object of type 'ref' has no len()
*/
%feature("python:slot", "tp_iter", functype="getiterfunc") ref_contig_iter::rapi___iter__;
%feature("python:slot", "tp_iternext", functype="iternextfunc") ref_contig_iter::next;
%extend ref_contig_iter {
  ref_contig_iter(rapi_contig* array, size_t len) {
    ref_contig_iter* iter = (ref_contig_iter*) rapi_malloc(sizeof(ref_contig_iter));
    if (!iter) return NULL;

    iter->next_item = array;
    iter->n_left = len;
    return iter;
  }

  ~ref_contig_iter() {
    free($self);
  }

  ref_contig_iter* rapi___iter__() { return $self; }

  rapi_contig* next() {
    if ($self->n_left > 0) {
      $self->n_left -= 1;
      return $self->next_item++;
    }
    else {
      PyErr_SetString(PyExc_StopIteration, "");
      return NULL;
    }
  }
};

/*
We're wrapping rapi_ref as a container of contigs.  In other words, the rapi_contig
elements are accessed by indexing into the rapi_ref.  Though I'd prefer it, I haven't
figured out how to expose the contigs attribute as an array of structures, supporting
usage like:

    > ref = rapi.ref('some path')
    > c1 = ref.contigs[0]
*/
%feature("python:slot", "sq_length", functype="lenfunc") rapi_ref::rapi___len__;
%feature("python:slot", "mp_subscript", functype="binaryfunc") rapi_ref::rapi___getitem__;
typedef struct {
  char * path;
  int n_contigs;
} rapi_ref;


%exception rapi_ref::rapi_get_contig {
  $action;
  if (result == NULL) {
    // The called $action should have already set the exception type and message
    SWIG_fail;
  }
}

%exception rapi_ref::rapi___getitem__ {
  $action;
  if (result == NULL) {
    // The called $action should have already set the exception type and message
    SWIG_fail;
  }
}

%feature("python:slot", "tp_iter", functype="getiterfunc") rapi_ref::rapi___iter__;
%extend rapi_ref {
  rapi_ref(const char* reference_path) {
    rapi_ref* ref = (rapi_ref*) rapi_malloc(sizeof(rapi_ref));
    if (!ref) return NULL;

    int error = rapi_ref_load(reference_path, ref);
    if (error == RAPI_NO_ERROR)
      return ref;
    else {
      free(ref);
      const char* msg;
      if (error == RAPI_MEMORY_ERROR)
        msg = "Insufficient memory available to load the reference.";
      else if (error == RAPI_GENERIC_ERROR)
        msg = "Library failed to load reference. Check paths and reference format.";
      else
        msg = "";

      PyErr_SetString(rapi_py_error_type(error), msg);
      return NULL;
    }
  }

  void unload() {
    rapi_error_t error = rapi_ref_free($self);
    if (error != RAPI_NO_ERROR) {
      PERROR("Problem destroying reference (error code %d)\n", error);
    }
  }

  ~rapi_ref() {
    rapi_ref_unload($self);
  }

  size_t rapi___len__() { return $self->n_contigs; }

  /* XXX: I worry about memory management here.  We're returning a pointer to the rapi_ref's
    chunk of memory.  I don't think there's anything preventing the interpreter from deciding
    to free the underlying memory and making everything blow up.
  */
  rapi_contig* rapi_get_contig(int i) {
    if (i < 0 || i >= $self->n_contigs) {
      PyErr_SetString(PyExc_IndexError, "index out of bounds");
      return NULL;
    }
    else
      return $self->contigs + i;
  }

  rapi_contig* rapi___getitem__(int i) { return rapi_ref_rapi_get_contig($self, i); }

  ref_contig_iter* rapi___iter__() { return new_ref_contig_iter($self->contigs, $self->n_contigs); }
};

/****** rapi_batch and rapi_reads *******/

%feature("python:slot", "sq_length", functype="lenfunc") rapi_read::rapi___len__;
typedef struct {
  char * id;
  char * seq;
  char * qual;
  unsigned int length;
  rapi_alignment* alignments;
  uint8_t n_alignments;
} rapi_read;

%extend rapi_read {
  size_t rapi___len__() { return $self->length; }
};



/******************************************
 * start:api_batch_wrap
 ******************************************/

/*
 * We wrap the rapi_batch in a shell class to implement some higher level
 * functionality -- principally appending.  To append we need to keep track
 * of how many reads have already been inserted (not just of allocated capacity,
 * which is what the bare rapi_batch C structure provides.  To add this new
 * member variable we have use this strategy.
 *
 * Notice that this structure is name "read_batch" in the wrapper interface
 * (thanks to a %rename rule earlier in this file).
*/

%feature("python:slot", "sq_length", functype="lenfunc") rapi_batch_wrap::rapi___len__;
%inline %{
  typedef struct {
    rapi_batch* batch;
    size_t len; // number of reads inserted in batch (as opposed to the space reserved)
  } rapi_batch_wrap;
%}

%extend rapi_batch_wrap {

  rapi_batch_wrap(int n_reads_per_frag) {
    if (n_reads_per_frag <= 0) {
      PyErr_SetString(PyExc_ValueError, "number of reads per fragment must be greater than or equal to 0");
      return NULL;
    }

    rapi_batch_wrap* wrapper = (rapi_batch_wrap*) rapi_malloc(sizeof(rapi_batch_wrap));
    if (!wrapper) return NULL;

    wrapper->batch = (rapi_batch*) rapi_malloc(sizeof(rapi_batch));
    if (!wrapper->batch) {
      free(wrapper);
      return NULL;
    }

    wrapper->len = 0;

    int error = rapi_reads_alloc(wrapper->batch, n_reads_per_frag, 0); // zero-sized allocation to initialize

    if (error != RAPI_NO_ERROR) {
      free(wrapper->batch);
      free(wrapper);
      PyErr_SetString(rapi_py_error_type(error), "Error allocating space for reads");
      return NULL;
    }
    else
      return wrapper;
  }

  ~rapi_batch_wrap() {
    int error = rapi_reads_free($self->batch);
    free($self);
    if (error != RAPI_NO_ERROR) {
      PERROR("Problem destroying read batch (error code %d)\n", error);
      // TODO: should we raise exceptions in case of errors when freeing/destroying?
    }
  }

  /** Number of reads per fragment */
  int n_reads_per_frag() { return $self->batch->n_reads_frag; }

  /** Number of reads inserted in batch (as opposed to the space reserved).
   *  This is actually index + 1 of the "forward-most" read to have been inserted.
   */
  int rapi___len__() { return $self->len; }

  /** Number of reads for which we have allocated memory. */
  int capacity() { return $self->batch->n_frags * $self->batch->n_reads_frag; }

  rapi_read* get_read(int n_fragment, int n_read)
  {
    //return rapi_get_read($self->batch, n_fragment, n_read);
    rapi_read* r = rapi_get_read($self->batch, n_fragment, n_read);
    PDEBUG("number of alignments in this read: %d\n", r->n_alignments);
    return r;
  }
}

%extend rapi_batch_wrap {

  rapi_error_t reserve(int n_reads) {
    if (n_reads < 0) {
      PyErr_SetString(rapi_py_error_type(RAPI_PARAM_ERROR), "number of reads to reserve must be >= 0");
      return RAPI_PARAM_ERROR;
    }

    int n_fragments = n_reads / $self->batch->n_reads_frag;
    // If the reads don't fit completely in n_fragments, add one more
    if (n_reads % $self->batch->n_reads_frag != 0)
      n_fragments +=  1;

    rapi_error_t error = rapi_reads_reserve($self->batch, n_fragments);

    if (error != RAPI_NO_ERROR) {
      PyErr_SetString(rapi_py_error_type(error), "Failed to reserve space");
    }
    /*else {
      for (int i = 0; i < $self->len; ++i) {
        PDEBUG("i: %d; n_alignments: %d\n", i, $self->batch->reads[i].n_alignments);
      }
    }*/
    return error;
  }

  rapi_error_t append(const char* id, const char* seq, const char* qual, int q_offset)
  {
    int fragment_num = $self->len / $self->batch->n_reads_frag;
    int read_num = $self->len % $self->batch->n_reads_frag;
    rapi_error_t error = RAPI_NO_ERROR;

    size_t read_capacity = rapi_batch_wrap_capacity($self);
    if (read_capacity < $self->len + 1) {
      // double the space
      size_t new_capacity = read_capacity > 0 ? read_capacity * 2 : 2;
      error = rapi_batch_wrap_reserve($self, new_capacity);
      if (error != RAPI_NO_ERROR) {
        return error;
      }
    }
    error = rapi_set_read($self->batch, fragment_num, read_num, id, seq, qual, q_offset);
    if (error != RAPI_NO_ERROR) {
      PyErr_SetString(rapi_py_error_type(error), "Error inserting read.");
    }
    else
      ++$self->len;

    return error;
  }

  int set_read(int n_frag, int n_read, const char* id, const char* seq, const char* qual, int q_offset)
  {
    if (n_frag < 0 || n_read < 0) {
      PyErr_SetString(rapi_py_error_type(RAPI_PARAM_ERROR), "read and fragment indices cannot be negative");
      return RAPI_PARAM_ERROR;
    }

    int error = rapi_set_read($self->batch, n_frag, n_read, id, seq, qual, q_offset);
    if (error != RAPI_NO_ERROR) {
      PyErr_SetString(rapi_py_error_type(error), "Error setting read data");
    } else {
      // If the user set a read that is beyond the current batch length, reset the
      // batch length to the new limit.
      int index = n_frag * $self->batch->n_reads_frag + n_read;
      if ($self->len < index)
        $self->len = index + 1;
    }
    return error;
  }
}

/******************************************
 * end:rapi_batch_wrap
 ******************************************/



%mutable;


// vim: set et sw=2 ts=2
