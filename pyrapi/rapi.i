

/* SWIG interface for rapi.h */

/******* Language-independent exceptions *******/
//%include "exception.i"
//SWIG_exception(SWIG_MemoryError, "Not enough memory");
// LP: I can't get SWIG_exception to work. It generates a "goto fail;"
// statement without defining the "fail" label, thus resulting in a
// compiler error.  I've tried this with SWIG 2.0.11 and 3.0.2.
// This leaves me writing one interface per target language.

/******* swig -builtin option *******
PyExc_SetString wouldn't raise a Python exception until I
started calling swig with the "-builtin" option.

*************************************/


%module rapi
%{
/* Includes the header in the wrapper code */
#include "rapi.h"
%}

%header %{

#include <stddef.h>
#include <stdio.h>

#define PDEBUG(...) { fprintf(stderr, "%s(%d): ", __FILE__, __LINE__); fprintf(stderr, __VA_ARGS__); }

/* Helper to convert RAPI errors to Python errors */
PyObject* rapi_py_error_type(int rapi_code) {
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

/*
Remove rapi_ prefix from function names.
Make sure this general rename comes before the ignore clauses below
or they won't have effect.
*/
%rename("%(strip:[rapi_])s") ""; // e.g., rapi_load_ref -> load_ref

// ignore functions from klib
%rename("$ignore", regextarget=1) "^k";
%rename("$ignore", regextarget=1) "^K";

// also ignore anything that starts with an underscore
%rename("$ignore", regextarget=1) "^_";

%ignore "rapi_init_kstr";

// opts are initialized when they're created.
%ignore "rapi_init_opts";

%include "kvec.h";
%include "kstring.h";
%include "rapi.h";

%extend rapi_opts {
  rapi_opts() {
    rapi_opts* opts = (rapi_opts*) rapi_malloc(sizeof(rapi_opts));
    if (!opts) return NULL;

    int error = rapi_init_opts(opts);
    if (error == RAPI_NO_ERROR)
      return opts;
    else {
      free(opts);
      PyErr_SetString(rapi_py_error_type(error), "Not enough memory to initialize rapi_opts");
      return NULL;
    }
  }

  ~rapi_opts() {
    int error = rapi_free_opts($self);
    if (error != RAPI_NO_ERROR) {
      PDEBUG("Problem destroying opts (error code %d)\n", error);
      // TODO: should we raise exceptions in case of errors when freeing/destroying?
    }
  }
};

// vim: set et sw=2 ts=2
