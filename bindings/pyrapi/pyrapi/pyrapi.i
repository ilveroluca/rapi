

/* SWIG interface for rapi.h */

%include "exception.i"


/******* swig -builtin option *******
PyErr_SetString wouldn't raise a Python exception until I
started calling swig with the "-builtin" option.

*************************************/


%module rapi

%include "../rapi_swig_helpers.i"

%header %{

/* includes injected into the C wrapper code.  */
#include <stddef.h>
#include <stdio.h>
#include <rapi.h>
#include <rapi_utils.h>

/* Then we define some helpers convert RAPI errors to swig errors */
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

/* malloc wrapper that sets a SWIG error in case of failure. */
void* rapi_malloc(size_t nbytes) {
  void* result = malloc(nbytes);
  if (!result)
    SWIG_Error(SWIG_MemoryError, "Failed to allocate memory");

  return result;
}

%}

/**** Now we begin wrapping things with SWIG */

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

/*** typemaps rapi_ssize_t <--> Python integer */
%typemap(out) rapi_ssize_t {
    $result = PyLong_FromLongLong($1);
    if (!$result) {
        SWIG_fail;
    }
}

%typemap(in) rapi_ssize_t {
    if (PyInt_Check($input)) {
        $1 = (rapi_ssize_t)PyInt_AsLong($input);
    }
    else if (PyLong_Check($input)) {
        $1 = (rapi_ssize_t)PyLong_AsLongLong($input);
    }
    else {
        SWIG_exception_fail(SWIG_TypeError, "Argument must be an integer");
    }

    if ($1 == -1 && PyErr_Occurred()) {
        // There was some other problem with the conversion
        SWIG_fail;
    }
}

/*
We take an analogous strategy to define a type rapi_bool that we only use
for this interface.
*/

%inline %{
typedef int rapi_bool;
%}

%typemap(out) rapi_bool {
    if ($1)
        $result = Py_True;
    else
        $result = Py_False;
    Py_INCREF($result);
}

/****
Define renaming and ignore rules.  The rules need to be defined
before the objects to which they are applied are declared to SWIG.
****/

// We don't provide direct access to the read_batch struct
%rename("$ignore") rapi_batch;
// Instead, we'll define a wrapper for it.  It's implemented
// as the structure rapi_batch_wrap within this wrapper and
// exposed to the user as rapi.read_batch
%rename(read_batch) rapi_batch_wrap;

// We'll also rename rapi_aligner_state to a simple "aligner"
%rename(aligner) rapi_aligner_state;

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


/*****
Tell swig about stdint and kvec
*****/
%include "stdint.i"; // needed to tell swig about types such as uint32_t, uint8_t, etc.

// things common to all bindings
%include "../rapi_common.i"

/***** A macro for an array iterator type*******/
%define rapi_array_iter(type)
%inline %{
typedef struct {
    const type* next_item;
    size_t n_left;
} type##_iter;
%}


/*
Without these %feature statements the "special" methods, such as __len__,
__getitem__, __iter__, etc. will be defined but they won't work (at least
when calling swig with -builtin).  You'll get errors like:

----> e.g., len(r)
TypeError: object of type 'ref' has no len()
*/
%feature("python:slot", "tp_iter", functype="getiterfunc") type##_iter::rapi___iter__;
%feature("python:slot", "tp_iternext", functype="iternextfunc") type##_iter::next;
%extend type##_iter {
    type##_iter(const type* array, size_t len) {
        if (array == NULL) {
            SWIG_Error(SWIG_ValueError, "array cannot be None");
            return NULL;
        }

        type##_iter* iter = (type##_iter*) rapi_malloc(sizeof(type##_iter));
        if (!iter) return NULL;
        iter->next_item = array;
        iter->n_left = len;
        return iter;
    }

    ~type##_iter(void) { free($self); }

    type##_iter* rapi___iter__(void) { return $self; }

    const type* next(void) {
        const type* retval;
        if ($self->n_left) {
            retval = $self->next_item++;
            $self->n_left--;
        }
        else {
          PyErr_SetString(PyExc_StopIteration, "");
          retval = NULL;
        }

        return retval;
    }
};
%enddef

/************ begin wrapping functions and structures **************/

/*
The char* returned by the following functions are wrapped automatically by
SWIG -- the wrapper doesn't try to free the strings since they are "const".
*/
const char* rapi_aligner_name(void);
const char* rapi_aligner_version(void);
const char* rapi_plugin_version(void);

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

  rapi_param(void) {
    rapi_param* param = (rapi_param*) rapi_malloc(sizeof(rapi_param));
    if (!param) return NULL;
    rapi_param_init(param);
    return param;
  }

  ~rapi_param(void) {
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
  int n_threads;
  rapi_bool share_ref_mem;

  /* Mismatch / Gap_Opens / Quality Trims --> Generalize ? */

  // TODO: how to wrap this thing?
  //kvec_t(rapi_param) parameters;
} rapi_opts;


%extend rapi_opts {
  rapi_opts(void) {
    rapi_opts* opts = (rapi_opts*) rapi_malloc(sizeof(rapi_opts));
    if (!opts) return NULL;

    rapi_error_t error = rapi_opts_init(opts);
    if (error == RAPI_NO_ERROR)
      return opts;
    else {
      free(opts);
      SWIG_Error(rapi_swig_error_type(error), "Not enough memory to initialize rapi_opts");
      return NULL;
    }
  }

  ~rapi_opts(void) {
    rapi_error_t error = rapi_opts_free($self);
    if (error != RAPI_NO_ERROR) {
      PERROR("Problem destroying opts (error code %d)\n", error);
      // TODO: should we raise exceptions in case of errors when freeing/destroying?
    }
  }
};

rapi_error_t rapi_init(const rapi_opts* opts);
rapi_error_t rapi_shutdown(void);


/****** IMMUTABLE *******/
// Everything from here down is read-only

%immutable;
/****** IMMUTABLE *******/

/***************************************
 ****** rapi_contig and rapi_ref *******
 ***************************************/

typedef struct {
  char * name;
  uint32_t len;
  char * assembly_identifier;
  char * species;
  char * uri;
  char * md5;
} rapi_contig;

/* An iterator object for the contig */
rapi_array_iter(rapi_contig);

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
    if (reference_path == NULL) {
        SWIG_Error(SWIG_TypeError, "Reference path cannot be None");
        return NULL;
    }

    rapi_ref* ref = (rapi_ref*) rapi_malloc(sizeof(rapi_ref));
    if (!ref) return NULL;

    rapi_error_t error = rapi_ref_load(reference_path, ref);
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

      SWIG_Error(rapi_swig_error_type(error), msg);
      return NULL;
    }
  }

  void unload(void) {
    rapi_error_t error = rapi_ref_free($self);
    if (error != RAPI_NO_ERROR) {
      PERROR("Problem destroying reference (error code %d)\n", error);
    }
  }

  ~rapi_ref(void) {
    rapi_ref_unload($self);
  }

  size_t rapi___len__(void) const { return $self->n_contigs; }

  /* XXX: I worry about memory management here.  We're returning a pointer to the rapi_ref's
    chunk of memory.  I don't think there's anything preventing the interpreter from deciding
    to free the underlying memory and making everything blow up.
  */
  const rapi_contig* rapi_get_contig(int i) const {
    if (i < 0 || i >= $self->n_contigs) {
      SWIG_Error(SWIG_IndexError, "index out of bounds");
      return NULL;
    }
    else
      return $self->contigs + i;
  }

  const rapi_contig* rapi___getitem__(int i) const { return rapi_ref_rapi_get_contig($self, i); }

  rapi_contig_iter* rapi___iter__(void) const { return new_rapi_contig_iter($self->contigs, $self->n_contigs); }
};

/***************************************
 ****** rapi_alignment           *******
 ***************************************/

// List of cigar ops.  We use the same strategy as for rapi_fragment
%typemap(out) rapi_cigar_ops {
    /**
    * Convert rapi_cigar_ops to tuples (char, len), where the operation
    * character corresponds to the one used in CIGAR notation:
    * 
    * M  0
    * I  1
    * D  2
    * S  3
    * H  4
    * N  5
    * P  6
    */
    // $1 is variable of tupel rapi_cigar_ops
    PyObject* ret_tuple = PyTuple_New($1.len);
    if (ret_tuple == NULL) {
        PERROR("Failed to allocate PyTuple. Raising exception\n");
        SWIG_exception_fail(SWIG_MemoryError, "Failed to allocate tuple");
    }

    // for each cigar operator in the array, map it to a 2-item tuple and append it to the return tuple
    for (int i = 0; i < $1.len; ++i) {
        rapi_cigar cig = $1.ops[i];
        char c = rapi_cigops_char[cig.op];
        long len = cig.len;

        PyObject* cigar_tuple = PyTuple_New(2);
        if (cigar_tuple == NULL) {
            Py_DECREF(ret_tuple);
            SWIG_exception_fail(SWIG_MemoryError, "Failed to allocate cigar tuple");
        }

        PyTuple_SET_ITEM(cigar_tuple, 0, PyString_FromStringAndSize(&c, 1));
        PyTuple_SET_ITEM(cigar_tuple, 1, PyInt_FromLong(len));
        PyTuple_SET_ITEM(ret_tuple, i, cigar_tuple);
    }

    $result = ret_tuple;
}

%{
typedef struct {
    rapi_cigar* ops;
    size_t len;
} rapi_cigar_ops;
%}

typedef struct {
} rapi_cigar_ops;

%{
/*
 * Wrap a tag value in a PyObject.
 *
 * Generates the appropriate type of Python object based on the tag value type.
 *
 * \return Sets exception and returns NULL in case of error.
 */
static PyObject* rapi_py_tag_value(const rapi_tag*const tag)
{
    if (tag == NULL) {
        SWIG_Error(SWIG_RuntimeError, "rapi_py_tag_value: Got NULL tag!");
        return NULL;
    }

    PyObject* retval = NULL;
    rapi_error_t error = RAPI_NO_ERROR;

    switch (tag->type) {
        case RAPI_VTYPE_CHAR: {
            char c;
            error = rapi_tag_get_char(tag, &c);
            if (error) {
                PERROR("rapi_tag_get_char returned error %s (%d)\n", rapi_error_name(error), error);
            }
            else
                retval = PyString_FromStringAndSize(&c, 1);
            break;
        }
        case RAPI_VTYPE_TEXT: {
            const kstring_t* s;
            error = rapi_tag_get_text(tag, &s);
            if (error) {
                PERROR("rapi_tag_get_text returned error %s (%d)\n", rapi_error_name(error), error);
            }
            else
                retval = PyString_FromStringAndSize(s->s, s->l);
            break;
        }
        case RAPI_VTYPE_INT: {
            long i;
            error = rapi_tag_get_long(tag, &i);
            if (error) {
                PERROR("rapi_tag_get_long returned error %s (%d)\n", rapi_error_name(error), error);
            }
            else
                retval = PyInt_FromLong(i);
            break;
        }
        case RAPI_VTYPE_REAL: {
            double d;
            error = rapi_tag_get_dbl(tag, &d);
            if (error) {
                PERROR("rapi_tag_get_dbl returned error %s (%d)\n", rapi_error_name(error), error);
            }
            else
                retval = PyFloat_FromDouble(d);
            break;
        }
        default:
            PERROR("Unrecognized tag type id %d\n", tag->type);
    }

    if (error != RAPI_NO_ERROR) {
        // this handler is for all the rapi_tag_get_* calls.  On the other hand, the Py functions
        // will set their own exception in case of error.
        SWIG_Error(rapi_swig_error_type(error), "rapi_py_tag_value: error extracting tag value");
        return NULL;
    }
    // else
    return retval;
}
%}

%typemap(out) rapi_tag_list {
    /**
     *  Map a rapi_tag_list to a dict( str -> value ).
     *  The key is always a string, but the value will be converted to the
     *  appropriate Python type.
     */
    PyObject* dict = PyDict_New();
    if (dict == NULL) {
        SWIG_exception_fail(SWIG_MemoryError, "Failed to allocate dict");
    }

    const char* error_msg = NULL;

    for (int i = 0; i < kv_size($1); ++i)
    {
        const rapi_tag*const pTag = &kv_A($1, i);
        PyObject* value = rapi_py_tag_value(pTag);
        if (value != NULL) {
            if (PyDict_SetItemString(dict, pTag->key, value) < 0)
                error_msg = "Error inserting tag into dict";
        }
        else
            error_msg = "Error converting tag to a python object";

        Py_XDECREF(value);
        if (error_msg)
            break;
    }

    if (error_msg) {
        Py_DECREF(dict);
        SWIG_exception_fail(SWIG_RuntimeError, error_msg);
    }
    else
        $result = dict;
};

typedef struct {
    rapi_contig* contig;
    unsigned long int pos; // 1-based
    uint8_t mapq;
    int score; // aligner-specific score

    uint8_t n_mismatches;
    uint8_t n_gap_opens;
    uint8_t n_gap_extensions;
} rapi_alignment;

%newobject rapi_alignment::get_cigar_string;
%extend rapi_alignment {

    // synthesize boolean attributes corresponding to bitfield values
    rapi_bool paired;
    rapi_bool prop_paired;
    rapi_bool mapped;
    rapi_bool reverse_strand;
    rapi_bool secondary_aln;

    rapi_cigar_ops get_cigar_ops(void) const {
        rapi_cigar_ops array;
        array.ops = $self->cigar_ops;
        array.len = $self->n_cigar_ops;
        return array;
    }

    char* get_cigar_string(void) const {
        kstring_t output = { 0, 0, NULL };
        rapi_put_cigar($self->n_cigar_ops, $self->cigar_ops, 0, &output);
        // return the string directly to Python who will be responsible for freeing it
        return output.s;
    }

    rapi_tag_list get_tags(void) const {
        return $self->tags;
    }

    int get_rlen(void) const {
      return rapi_get_rlen($self->n_cigar_ops, $self->cigar_ops);
    }
};


%{
rapi_bool rapi_alignment_paired_get(const rapi_alignment* aln) {
    return aln->paired != 0;
}

rapi_bool rapi_alignment_prop_paired_get(const rapi_alignment* aln) {
    return aln->prop_paired != 0;
}

rapi_bool rapi_alignment_mapped_get(const rapi_alignment* aln) {
    return aln->mapped != 0;
}

rapi_bool rapi_alignment_reverse_strand_get(const rapi_alignment* aln) {
    return aln->reverse_strand != 0;
}

rapi_bool rapi_alignment_secondary_aln_get(const rapi_alignment* aln) {
    return aln->secondary_aln != 0;
}
%}



/***************************************
 ****** rapi_batch and rapi_reads ******
 ***************************************/

/********************** typemap rapi_fragment -> Tuple ******************/
/*
  Map a rapi_fragment to a Python tuple containing the relevant reads.
  NOTE:
  We don't copy any data here.  We merely create SWIG wrappers for the
  rapi_reads returned by rapi_get_read(batch, f, r).  Therefore,
  the Python object will have 'thisown' set to 0 AND the user must keep
  a reference to the read batch around while he uses the data.
*/
%typemap(out) rapi_fragment {

  int fragment_size = $1.batch->n_reads_frag;
  PyObject* tuple = PyTuple_New(fragment_size);
  if (tuple == NULL) {
    PERROR("Failed to allocate PyTuple. Raising exception\n");
    SWIG_exception_fail(SWIG_MemoryError, "Failed to allocate tuple");
  }

  for (int i = 0; i < fragment_size; ++i) {
    const rapi_read* read = rapi_get_read($1.batch, $1.fragment_num, i);
    if (read == NULL) {
      PERROR("Failed to fetch read (rapi_get_read(batch, %lld, %d) returned NULL\n", $1.fragment_num, i);
      Py_DECREF(tuple);
      SWIG_exception_fail(SWIG_RuntimeError, "Error retrieving read");
    }
    // SWIG_NewPointerObj creates a new Python wrapper for the rapi_read to which
    // we're pointing. By passing the flag 0 we're telling SWIG that the wrapper
    // doesn't own the rapi_read, thus it should not try to free it.
    PyObject* read_object = SWIG_NewPointerObj(SWIG_as_voidptr(read), SWIGTYPE_p_rapi_read, 0);
    // PyTuple_SET_ITEM steals the object reference so we don't need to DECREF.
    PyTuple_SET_ITEM(tuple, i, read_object);
  }
  $result = tuple;
}

// define rapi_alignment_iter: iterator for array of rapi_alignment
rapi_array_iter(rapi_alignment);

// Raise IndexError when rapi_read.get_aln returns NULL
%exception rapi_read::get_aln {
    $action
    if (result == NULL) {
        SWIG_exception(SWIG_IndexError, "");
    }
}

%feature("python:slot", "sq_length", functype="lenfunc") rapi_read::rapi___len__;
typedef struct {
    char * id;
    char * seq;
    char * qual;
    unsigned int length;
    uint8_t n_alignments;
} rapi_read;

%extend rapi_read {
    size_t rapi___len__(void) const { return $self->length; }

    const rapi_alignment* get_aln(int index) const {
        if (index >= 0 && index < $self->n_alignments)
            return $self->alignments + index;
        else
            return NULL; // exception raise in %exception block
    }

    rapi_alignment_iter* iter_aln(void) const {
        return new_rapi_alignment_iter($self->alignments, $self->n_alignments);
    }

    // synthesize some high-level boolean attributes.  These
    // look at the "main" alignment (first one, if any)
    rapi_bool prop_paired;
    rapi_bool mapped;
    rapi_bool reverse_strand;
    int score;
    uint8_t mapq;
};

%{
rapi_bool rapi_read_prop_paired_get(const rapi_read* read) {
    return read->n_alignments > 0 && rapi_alignment_prop_paired_get(read->alignments);
}

rapi_bool rapi_read_mapped_get(const rapi_read* read) {
    return read->n_alignments > 0 && rapi_alignment_mapped_get(read->alignments);
}

rapi_bool rapi_read_reverse_strand_get(const rapi_read* read) {
    return read->n_alignments > 0 && rapi_alignment_reverse_strand_get(read->alignments);
}

int rapi_read_score_get(const rapi_read* read) {
    if (read->n_alignments <= 0) {
        return 0;
    }
    else {
        return read->alignments->score;
    }
}

uint8_t rapi_read_mapq_get(const rapi_read* read) {
    if (read->n_alignments <= 0) {
        return 0;
    }
    else {
        return read->alignments->mapq;
    }
}

%}


/*
 * We wrap the rapi_batch in a shell class to implement some higher level
 * functionality -- principally appending.  To append we need to keep track
 * of how many reads have already been inserted (not just of allocated capacity,
 * which is what the bare rapi_batch C structure provides.  To add this new
 * member variable we have use this strategy.
 *
 * NOTE: this structure is called *read_batch* in the wrapper interface
 * (thanks to a %rename rule earlier in this file).
*/

// Forward declarations
%{
struct rapi_batch_wrap;
struct read_batch_iter;
SWIGINTERN struct read_batch_iter* new_read_batch_iter(const struct rapi_batch_wrap* batch);
%}

%feature("python:slot", "sq_length", functype="lenfunc") rapi_batch_wrap::rapi___len__;
%feature("python:slot", "tp_iter", functype="getiterfunc") rapi_batch_wrap::rapi___iter__;
%{ // this declaration is inserted in the C code
typedef struct rapi_batch_wrap {
  rapi_batch* batch;
  rapi_ssize_t len; // number of reads inserted in batch (as opposed to the space reserved)
} rapi_batch_wrap;
%}

/*
 * Implement two synthesized attributes for this object.  These functions
 * are up here so that in the generated code the appear before their first
 * invocation.
 */
%{
int rapi_batch_wrap_n_reads_per_frag_get(const rapi_batch_wrap* self) {
    return self->batch->n_reads_frag;
}

rapi_ssize_t rapi_batch_wrap_n_fragments_get(const rapi_batch_wrap* self) {
    return self->len / self->batch->n_reads_frag;
}

rapi_ssize_t rapi_batch_wrap_capacity_get(const rapi_batch_wrap* wrap) {
  return rapi_batch_read_capacity(wrap->batch);
}

%}

// This one to the SWIG interpreter.
// We don't expose any of the struct members through SWIG.
typedef struct rapi_batch_wrap {
} rapi_batch_wrap;

%{
typedef struct rapi_fragment {
  const rapi_batch* batch;
  rapi_ssize_t fragment_num;
} rapi_fragment;
%}
typedef struct rapi_fragment {
} rapi_fragment;

/* Iterator over fragments in a read batch */
// don't expose (%ignore) the members of the read_batch_iter struct
%ignore read_batch_iter::wrapper;
%ignore read_batch_iter::next_fragment;
%feature("python:slot", "tp_iter", functype="getiterfunc") read_batch_iter::__iter__;
%feature("python:slot", "tp_iternext", functype="iternextfunc") read_batch_iter::next;
%{
typedef struct read_batch_iter {
  const rapi_batch_wrap* wrapper;
  size_t next_fragment;
} read_batch_iter;
%}
typedef struct read_batch_iter {
} read_batch_iter;

%exception rapi_batch_wrap::get_read {
  $action
  if (result == NULL) {
    SWIG_exception_fail(SWIG_IndexError, "co-ordinates out of bounds");
  }
}

%extend rapi_batch_wrap {

  /**
   * Creates a new read_batch for fragments composed of `n_reads_per_frag` reads.
   * The function doesn't pre-allocate any space for reads, so either use `append`
   * to insert reads or call `reserve` before calling `set_read`.
   */
  rapi_batch_wrap(int n_reads_per_frag) {
    if (n_reads_per_frag <= 0) {
      SWIG_Error(SWIG_ValueError, "number of reads per fragment must be greater than or equal to 0");
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

    rapi_error_t error = rapi_reads_alloc(wrapper->batch, n_reads_per_frag, 0); // zero-sized allocation to initialize

    if (error != RAPI_NO_ERROR) {
      free(wrapper->batch);
      free(wrapper);
      SWIG_Error(rapi_swig_error_type(error), "Error allocating space for reads");
      return NULL;
    }
    else
      return wrapper;
  }

  ~rapi_batch_wrap(void) {
    rapi_error_t error = rapi_reads_free($self->batch);
    free($self);
    if (error != RAPI_NO_ERROR) {
      PERROR("Problem destroying read batch (error code %d)\n", error);
      // TODO: should we raise exceptions in case of errors when freeing/destroying?
    }
  }

  /** Number of reads per fragment */
  const int n_reads_per_frag;

  /** Number of complete fragments inserted */
  const rapi_ssize_t n_fragments;

  /** Number of reads for which we have allocated memory. */
  const rapi_ssize_t capacity;

  /** Number of reads inserted in batch (as opposed to the space reserved).
   *  This is actually index + 1 of the "forward-most" read to have been inserted.
   */
  rapi_ssize_t rapi___len__(void) const { return $self->len; }

  const rapi_read* get_read(rapi_ssize_t n_fragment, int n_read) const
  {
    // Since the underlying code merely checks whether we're indexing
    // allocated space, we precede it with an additional check whether we're
    // accessing space where reads have been inserted.
    //
    // * In both cases, if the returned value is NULL an IndexError is raised
    // in the %exception block
    if (n_fragment * $self->batch->n_reads_frag + n_read >= $self->len) {
        return NULL;
    }

    return rapi_get_read($self->batch, n_fragment, n_read);
  }

  rapi_error_t reserve(rapi_ssize_t n_reads) {
    if (n_reads < 0) {
        PERROR("n_reads must be >= 0");
        return RAPI_PARAM_ERROR;
    }

    rapi_ssize_t n_fragments = n_reads / $self->batch->n_reads_frag;
    // If the reads don't fit completely in n_fragments, add one more
    if (n_reads % $self->batch->n_reads_frag != 0)
      n_fragments +=  1;

    rapi_error_t error = rapi_reads_reserve($self->batch, n_fragments);
    return error;
  }

  rapi_error_t append(const char* id, const char* seq, const char* qual, int q_offset)
  {
    rapi_error_t error = RAPI_NO_ERROR;

    // if id or seq are NULL set them to the empty string and pass them down to the plugin.
    if (!id) id = "";
    if (!seq) seq = "";

    rapi_ssize_t fragment_num = $self->len / $self->batch->n_reads_frag;
    int read_num = $self->len % $self->batch->n_reads_frag;

    rapi_ssize_t read_capacity = rapi_batch_wrap_capacity_get($self);
    if (read_capacity < $self->len + 1) {
      // double the space
      rapi_ssize_t new_capacity = read_capacity > 0 ? read_capacity * 2 : 2;
      error = rapi_batch_wrap_reserve($self, new_capacity);
      if (error != RAPI_NO_ERROR) {
        return error;
      }
    }
    error = rapi_set_read($self->batch, fragment_num, read_num, id, seq, qual, q_offset);
    if (error != RAPI_NO_ERROR) {
      SWIG_Error(rapi_swig_error_type(error), "Error inserting read.");
    }
    else
      ++$self->len;

    return error;
  }

  rapi_error_t clear(void) {
    rapi_error_t error = rapi_reads_clear($self->batch);
    if (error == RAPI_NO_ERROR)
      $self->len = 0;
    return error;
  }

  rapi_error_t set_read(rapi_ssize_t n_frag, int n_read, const char* id, const char* seq, const char* qual, int q_offset)
  {
    // if id or seq are NULL set them to the empty string and pass them down to the plugin.
    if (!id) id = "";
    if (!seq) seq = "";

    rapi_error_t error = rapi_set_read($self->batch, n_frag, n_read, id, seq, qual, q_offset);
    if (error != RAPI_NO_ERROR) {
      SWIG_Error(rapi_swig_error_type(error), "Error setting read data");
    }
    else {
      // If the user set a read that is beyond the current batch length, reset the
      // batch length to the new limit.
      rapi_ssize_t index = n_frag * $self->batch->n_reads_frag + n_read;
      if ($self->len <= index)
        $self->len = index + 1;
    }
    return error;
  }

  // need the 'struct' keyword in the function prototype because we haven't
  // defined the read_batch_iter struct yet.
  struct read_batch_iter* rapi___iter__(void) const {
    return new_read_batch_iter($self);
  }
}

%exception read_batch_iter::next {
  $action
  if (result.batch == NULL)
    SWIG_fail; // exception already set by call
}

%extend read_batch_iter {
  read_batch_iter(const rapi_batch_wrap * batch) {
    if (batch == NULL) {
        SWIG_Error(SWIG_TypeError, "batch cannot be None");
        return NULL;
    }

    read_batch_iter* iter = (read_batch_iter*) rapi_malloc(sizeof(read_batch_iter));
    if (iter == NULL) return NULL;

    iter->wrapper = batch;
    iter->next_fragment = 0;
    return iter;
  }

  ~read_batch_iter(void) {
    free($self);
  }

  read_batch_iter* __iter__(void) { return $self; }

  rapi_fragment next() {
    rapi_fragment fragment;

    if ($self->next_fragment < rapi_batch_wrap_n_fragments_get($self->wrapper)) {
      fragment.batch = $self->wrapper->batch;
      fragment.fragment_num = $self->next_fragment;
      $self->next_fragment += 1;
    }
    else {
      PyErr_SetString(PyExc_StopIteration, "");
      fragment.batch = NULL;
    }
    return fragment;
  }
}

/***************************************
 ****** rapi_aligner             *******
 ***************************************/

%{ // forward declaration of opaque structure (in C-code)
struct rapi_aligner_state;
%}

// declare the structure to SWIG as an empty struct
typedef struct {
} rapi_aligner_state;

// attach methods to it
%extend rapi_aligner_state {
  rapi_aligner_state(const rapi_opts* opts) {
    struct rapi_aligner_state* pState;
    rapi_error_t error = rapi_aligner_state_init(&pState, opts);
    if (error != RAPI_NO_ERROR) {
      SWIG_Error(rapi_swig_error_type(error), "Error initializing aligner");
      return NULL;
    }
    return pState;
  }

  ~rapi_aligner_state(void) {
    rapi_error_t error = rapi_aligner_state_free($self);
    if (error != RAPI_NO_ERROR)
      PERROR("Problem destroying aligner state object (error code %d)\n", error);
  }

  rapi_error_t align_reads(const rapi_ref* ref, rapi_batch_wrap* batch) {
    if (NULL == ref || NULL == batch) {
      PERROR("ref and batch arguments must not be NULL\n");
      return RAPI_PARAM_ERROR;
    }

    if (batch->len % batch->batch->n_reads_frag != 0) {
      PERROR("Incomplete fragment in batch! Number of reads appended (%lld) is not a multiple of the number of reads per fragment (%d)\n",
        batch->len, batch->batch->n_reads_frag);
      return RAPI_GENERIC_ERROR;
    }

    rapi_ssize_t start_fragment = 0;
    rapi_ssize_t end_fragment = batch->len / batch->batch->n_reads_frag;
    return rapi_align_reads(ref, batch->batch, start_fragment, end_fragment, $self);
  }
}


/***************************************
 ****** other stuff              *******
 ***************************************/
%rename(rev_comp) rapi_rev_comp_wrapper;

// The next two typemaps configure the OutValue for the rapi_rev_comp_wrapper
// function.
%typemap(argout) PyObject** OutValue {
  $result = *$1;
}
// This one tells SWIG that the OutValue argument isn't used for input.
// It also allocates a temporary local variable to hold the output value,
// and sets $1 to point to it.
%typemap(in, numinputs=0) PyObject** OutValue(PyObject* temp) {
  $1 = &temp;
}

rapi_error_t rapi_rev_comp_wrapper(PyObject* seq, PyObject** OutValue);

%{
rapi_error_t rapi_rev_comp_wrapper(PyObject* seq, PyObject** outRevComp) {
  if (!PyString_Check(seq)) {
    return RAPI_TYPE_ERROR;
  }

  Py_ssize_t seq_len = PyString_GET_SIZE(seq);
  PyObject* retval = PyString_FromStringAndSize(PyString_AS_STRING(seq), seq_len);
  if (!retval) {
    return RAPI_GENERIC_ERROR;
  }

  char* tmp_str = PyString_AS_STRING(retval);
  rapi_error_t error = rapi_rev_comp(tmp_str, seq_len);
  if (error != RAPI_NO_ERROR) {
    Py_DECREF(retval);
    retval = NULL;
  }

  *outRevComp = retval;
  return error;
}
%}

/***
 * Swig's default output typemap to wrap char* strings converts NULL
 * into a None (in Python).  The following functions return NULL to indicate
 * an error so we need to handle them differently.
 *
 * NOTE:  keep this stuff at the end of this file, unless you plan on restoring
 * the default output char* typemap.
 */
%typemap(out) char* {

  if ($1) {
    // convert char* to a string in the target language.
    // Unfortunately with Python this requires making a copy of the string.
    $result = SWIG_FromCharPtr((const char *)$1);
    free((char*)$1);
  }
  else {
    $result = NULL;
  }
}

%inline %{

/**
 * Input: some python sequence of rapi reads.
 */
char* format_sam(PyObject* pyobj) {
  char* retval = NULL;
  const rapi_read** read_ptrs = NULL;

  if (!PySequence_Check(pyobj)) {
    SWIG_Error(SWIG_TypeError, "Expected a sequence type");
    return NULL;
  }

  Py_ssize_t len = PySequence_Length(pyobj);
  if (len < 0) {
    SWIG_Error(SWIG_RuntimeError, "Failed to get length of sequence");
    return NULL;
  }

  read_ptrs = calloc(len, sizeof(rapi_read*));
  if (NULL == read_ptrs) {
    SWIG_Error(SWIG_MemoryError, "Failed to allocate memory");
    goto clean;
  }

  for (Py_ssize_t i = 0; i < len; ++i) {
    PyObject* swig_read = PySequence_Fast_GET_ITEM(pyobj, i); // borrowed reference
    // Get the rapi read from the python sequence and write a pointer to it into our `read_ptrs` array
    if (!swig_read
      || !SwigPyObject_Check(swig_read)
      || !SWIG_IsOK(SWIG_ConvertPtr(swig_read, (void**)&(read_ptrs[i]), SWIGTYPE_p_rapi_read, 0))) {
      SWIG_Error(SWIG_TypeError, "Expecting a sequence a RAPI read_ptrs");
      goto clean;
    }
  }

  // now we can finally call our formatting function
  kstring_t str = { 0, 0, NULL };
  rapi_error_t error = rapi_format_sam(read_ptrs, len, &str);

  if (error == RAPI_NO_ERROR)
    retval = str.s; // Python must free this string
  else {
    free(str.s);
    SWIG_Error(rapi_swig_error_type(error), "Error formatting SAM");
    // leave retval as NULL
  }

clean:
  free((rapi_read*)read_ptrs);
  return retval;
}


char* format_sam_from_batch(const rapi_batch_wrap* wrapper, rapi_ssize_t n_frag) {
  if (NULL == wrapper) {
    SWIG_Error(SWIG_TypeError, "wrapper argument cannot be None");
    return NULL;
  }

  if (n_frag < 0 || n_frag >= wrapper->len) {
    SWIG_Error(SWIG_IndexError, "n_frag argument out of bounds\n");
    return NULL;
  }

  kstring_t str = { 0, 0, NULL };
  rapi_error_t error = rapi_format_sam_b(wrapper->batch, n_frag, &str);
  if (error == RAPI_NO_ERROR)
    return str.s; // Python must free this string
  else {
    free(str.s);
    SWIG_Error(rapi_swig_error_type(error), "Error formatting SAM");
    return NULL;
  }
}

char* format_sam_hdr(const rapi_ref* ref)
{
  if (NULL == ref) {
    SWIG_Error(SWIG_TypeError, "ref argument cannot be None");
    return NULL;
  }

  kstring_t str = { 0, 0, NULL };
  rapi_error_t error = rapi_format_sam_hdr(ref, &str);
  if (error == RAPI_NO_ERROR)
    return str.s; // Python must free this string
  else {
    free(str.s);
    SWIG_Error(rapi_swig_error_type(error), "Error formatting SAM header");
    return NULL;
  }
}
%}

long rapi_get_insert_size(const rapi_alignment* read, const rapi_alignment* mate);

// vim: set et sw=2 ts=2
