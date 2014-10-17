


#define PDEBUG(...) { fprintf(stderr, "%s(%d) DEBUG: ", __FILE__, __LINE__); fprintf(stderr, __VA_ARGS__); }
#define PERROR(...) { fprintf(stderr, "%s(%d) ERROR: ", __FILE__, __LINE__); fprintf(stderr, __VA_ARGS__); }
