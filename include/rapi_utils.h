
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
******************************************************************************/

#ifndef __RAPI_UTILS_H__
#define __RAPI_UTILS_H__


#define PDEBUG(...) { fprintf(stderr, "%s(%d) DEBUG: ", __FILE__, __LINE__); fprintf(stderr, __VA_ARGS__); }
#define PERROR(...) { fprintf(stderr, "%s(%d) ERROR: ", __FILE__, __LINE__); fprintf(stderr, __VA_ARGS__); }

/* Utilities section */

long rapi_get_insert_size(const rapi_alignment* read, const rapi_alignment* mate);

/** The length of the aligned read in terms of reference bases. */
int rapi_get_rlen(int n_cigar, const rapi_cigar* cigar_ops);

/**
 * Write alignment operations in `ops` to a cigar string in `output`.
 *
 * \param ops Array of rapi_cigar.
 * \param n_opts Length of `ops`.
 * \param force_hard_clip If true, soft clips are converted to hard clips.
 * \param output String into which the output is appended.
 *
 * \note `output` must be initialized prior to calling this function.  The new
 * cigar string will be appended.
 */
void rapi_put_cigar(int n_ops, const rapi_cigar* ops, int force_hard_clip, kstring_t* output);


/******* SAM output *******/

/**
 * Format SAM for all reads in the given fragment.
 *
 * \param reads: pointer to list of read pointers; reads must be ordered first to last
 *
 * \param n_reads: number of reads in list
 *
 * \param output An initialized kstring_t to which the SAM will be appended.
 */
rapi_error_t rapi_format_sam(const rapi_read** reads, int n_reads, kstring_t* output);

/**
 * Format SAM for all reads in the indicated fragment.
 *
 * \param batch A read batch.
 *
 * \param n_frag The index of the fragment within the batch whose reads will be
 *               printed in the SAM text.  SAM will be produced for all the reads
 *               within the fragment.
 *
 * \param output An initialized kstring_t to which the SAM will be appended.
 */
rapi_error_t rapi_format_sam_b(const rapi_batch* batch, rapi_ssize_t n_frag, kstring_t* output);

/**
 * Format the SAM header for the given reference.  The header will also contain
 * a @PG tag identifying the RAPI-interfaced aligner being used.
 *
 * \param output An initialized kstring_t to which the output will be appended.
 */
rapi_error_t rapi_format_sam_hdr(const rapi_ref* ref, kstring_t* output);



/**
 * Compute the reverse complement of a sequence, in place.
 *
 * \param seq Array of uppercase characters in set {A, C, G, N, T }.
 * \param len Lengthof `seq`.
 *
 * \note In case of error, the sequence `seq` may be damaged (i.e., partially reverse-complemented.
 */
rapi_error_t rapi_rev_comp(char* seq, int len);

#endif
