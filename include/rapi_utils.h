
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
