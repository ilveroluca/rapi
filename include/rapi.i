

/* SWIG interface for rapi.h */

%module rapi
%{
/* Includes the header in the wrapper code */
#include "rapi.h"
%}

// ignore functions from klib
%rename("$ignore", regextarget=1) "^k";
%rename("$ignore", regextarget=1) "^K";
// also ignore anything that starts with an underscore
%rename("$ignore", regextarget=1) "^_";
%ignore "rapi_init_kstr";

%include "kvec.h"
%include "kstring.h"
%include "rapi.h"
