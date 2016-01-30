
#ifndef __RAPI_UTILS_H__
#define __RAPI_UTILS_H__


#define PDEBUG(...) { fprintf(stderr, "%s(%d) DEBUG: ", __FILE__, __LINE__); fprintf(stderr, __VA_ARGS__); }
#define PERROR(...) { fprintf(stderr, "%s(%d) ERROR: ", __FILE__, __LINE__); fprintf(stderr, __VA_ARGS__); }

rapi_error_t rapi_rev_comp(char* seq, int len);

#endif
