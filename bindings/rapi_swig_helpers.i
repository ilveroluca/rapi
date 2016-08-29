
/*************************************************
 * Import this file at the start of the language-specific
 * Swig interface definition.  Its contents need to go into
 * Swig's %header section, before other definitions are generated.
 */



// Default output typemap for char*
%typemap(out) char_ptr_default = char*;

/***
 * Swig's default output typemap to wrap char* strings converts
 * NULL into a None (in Python).  Some of our functions return NULL
 * to indicate an error so we need to handle them differently.  For
 * these we define an alternative output typemap for char*.
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

%typemap(out) char_ptr_null_error = char*;

// now restore default behaviour
%typemap(out) char* = char_ptr_default;



%header %{

/* includes injected into the C wrapper code.  */
#include <stddef.h>
#include <stdio.h>
#include <rapi.h>
#include <rapi_utils.h>

%}
