/*
 * rapi_bwa.c
 */

#include <rapi.h>
#include <rapi_utils.h>
#include <bwamem.h>
#include <kstring.h>
#include <kvec.h>
#include <utils.h>

#include <inttypes.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "bwa_header.h"

#define RAPI_BWA_PLUGIN_VERSION  "0.1.0-dev"

/* Internal configuration structure */
typedef struct {
	int mapq_min;
	int isize_min;
	int isize_max;
	int n_threads;
    mem_opt_t* bwa_opts;
} library_opts;

library_opts* _g_library_opts = NULL;

const char vtype_char[] = {
	'0',
	'A', // RAPI_VTYPE_CHAR       1
	'Z', // RAPI_VTYPE_TEXT       2
	'i', // RAPI_VTYPE_INT        3
	'f'  // RAPI_VTYPE_REAL       4
};

/**
 * The 'bwa_pg' string is statically allocated in some BWA file that we're not
 * linking, so we need to define it here.  I think it's the string that the
 * program uses to identify itself in the SAM header's @PG tag.
 */
const char*const bwa_pg = "rapi-bwa";

// Table to convert from base letter to an index in the range [0,4].
// Defined in bntseq.c
extern unsigned char nst_nt4_table[256];

/******** Utility functions *******/
void rapi_print_read(FILE* out, const rapi_read* read)
{
	fprintf(out, "read id: %s\n", read->id);
	fprintf(out, "read length: %d\n", read->length);
	fprintf(out, "read seq: %s\n", read->seq);
	fprintf(out, "read qual: %s\n", read->qual);
	fprintf(out, "read n_alignments: %u\n", read->n_alignments);
}

static void rapi_print_batch(FILE* out, const rapi_batch* batch)
{
	for (int f = 0; f < batch->n_frags; ++f) {
		for (int r = 0; r < batch->n_reads_frag; ++r) {
			const rapi_read*const read = rapi_get_read(batch, f, r);
			if (!read) {
				PERROR("rapi_get_read(batch, %d, %d) returned NULL\n", f, r);
				abort();
			}
			fprintf(out, "=================== (%d, %d) ===================\n ", f, r);
			rapi_print_read(out, read);
		}
	}
}

static void rapi_print_bwa_flag_string(FILE* out, const int flag)
{
	fprintf(out, "BWA flags: ");
	if (flag & MEM_F_PE)        fprintf(out, " MEM_F_PE");
	if (flag & MEM_F_NOPAIRING) fprintf(out, " MEM_F_NOPAIRING");
	if (flag & MEM_F_ALL)       fprintf(out, " MEM_F_ALL");
	if (flag & MEM_F_NO_MULTI)  fprintf(out, " MEM_F_NO_MULTI");
	if (flag & MEM_F_NO_RESCUE) fprintf(out, " MEM_F_NO_RESCUE ");
	if (flag & MEM_F_NO_EXACT)  fprintf(out, " MEM_F_NO_EXACT");
	fprintf(out, "\n");
}

/*
 * Writes a null-terminated string representation of the SAM flag to buf20.
 * Guaranteed not to write more than 20 bytes.
 *
 * \param buf20: ptr to character buffer at least 20 bytes long
 *
 * \returns the number of characters written, not including the NULL terminator.
 */
int rapi_flag_string(const int flag, char* buf20)
{
	const int max_length = 20 - 1; // 1 for null terminator
	char names[] = { 'p',    'P',    'u',    'U',    'r',    'R',    '1',    '2',    's',    'f',    'd',   's' };
	int values[] = { 0x0001, 0x0002, 0x0004, 0x0008, 0x0010, 0x0020, 0x0040, 0x0080, 0x0100, 0x0200, 0x0400, 0x0800 };
	int str_pos = 0;

	for (int i = 0; (i < sizeof(names) / sizeof(names[0])) && str_pos < max_length; ++i) {
		if (flag & values[i]) {
			buf20[str_pos] = names[i];
			str_pos += 1;
		}
	}
	buf20[str_pos] = '\0';
	return str_pos;
}

rapi_error_t rapi_format_tag(const rapi_tag* tag, kstring_t* str) {
	// in theory we should check the return values of all these kput functions
	// and ensure they're != EOF
	rapi_error_t error = RAPI_NO_ERROR;

	kputs(tag->key, str);
	kputc(':', str);
	kputc(vtype_char[tag->type], str);
	kputc(':', str);
	switch (tag->type) {
		case RAPI_VTYPE_CHAR: {
			char c;
			error = rapi_tag_get_char(tag, &c);
			if (error) return RAPI_TYPE_ERROR;
			kputc(c, str);
			break;
		 }
		case RAPI_VTYPE_TEXT: {
			const kstring_t* s;
			error = rapi_tag_get_text(tag, &s);
			if (error) return RAPI_TYPE_ERROR;
			kputsn(s->s, s->l, str);
			break;
		}
		case RAPI_VTYPE_INT: {
			long i;
			error = rapi_tag_get_long(tag, &i);
			if (error) return RAPI_TYPE_ERROR;
			kputl(i, str);
			break;
		}
		case RAPI_VTYPE_REAL: {
			double d;
			error = rapi_tag_get_dbl(tag, &d);
			if (error) return RAPI_TYPE_ERROR;
			ksprintf(str, "%f", d);
			break;
		}
		default:
			err_fatal(__func__, "Unrecognized tag type id %d\n", tag->type);
			abort();
	};
	return error;
}

/**
 * Produce SAM for `read`, using the alignment at index i_aln, or no alignment (as unmapped read) if i_aln < 0.
 */
static rapi_error_t _rapi_format_sam_aln(const rapi_read* read, int i_aln, const rapi_read* mate, int read_num, kstring_t* output)
{
	/**** code based on mem_aln2sam in BWA ***/

	if (NULL == read) {
		PERROR("_rapi_format_sam_aln: NULL read pointer\n");
		return RAPI_PARAM_ERROR;
	}

	if (read->n_alignments > 0 && i_aln >= read->n_alignments) {
		PERROR("_rapi_format_sam_aln: i_aln out of bounds\n");
		return RAPI_PARAM_ERROR;
	}

	rapi_alignment tmp_read, tmp_mate;

	if (i_aln < 0) { // select no alignment
		memset(&tmp_read, 0, sizeof(tmp_read));
	}
	else {
		tmp_read = read->alignments[i_aln];
	}

	if (mate && mate->n_alignments > 0) {
		tmp_mate = *mate->alignments;
	}
	else {
		memset(&tmp_mate, 0, sizeof(tmp_mate));
	}

	if (mate) {
		tmp_read.paired = 1;
		tmp_mate.paired = 1;
	}

	rapi_alignment* aln = &tmp_read;
	rapi_alignment* mate_aln = &tmp_mate;

	if (!aln->mapped && mate && mate_aln->mapped) { // copy mate position to read
		aln->contig         = mate_aln->contig;
		aln->pos            = mate_aln->pos;
		aln->reverse_strand = mate_aln->reverse_strand;
	}
	else if (aln->mapped && mate && !mate_aln->mapped) { // copy read alignment to mate
		mate_aln->contig         = aln->contig;
		mate_aln->pos            = aln->pos;
		mate_aln->reverse_strand = aln->reverse_strand;
	}

	int flag = 0;

	if (read_num == 1) flag |= 0x40;
	if (read_num == 2) flag |= 0x80;

	flag |= (mate && !mate_aln->mapped) ? 0x8 : 0; // is mate unmapped
	// for the 0x20 flag, we set it regardless of whether the mate is mapped.  If
	// the mate is not mapped but the read is, then the flag will have been copied
	// from the read a few lines above.  This replicates BWA's behaviour.
	flag |= (mate && mate_aln->reverse_strand) ? 0x20 : 0; // is mate on the reverse strand

	flag |= aln->paired ? 0x1 : 0; // is paired in sequencing
	flag |= aln->mapped ? 0 : 0x4; // is unmapped

	// As for the 0x20 flag, the same logic goes here for the 0x10 flag.
	// If the read isn't mapped aln->reverse_strand will have been copied from the mate.
	flag |= aln->reverse_strand ? 0x10 : 0; // is on the reverse strand

	if (aln->mapped)
	{
		flag |= aln->prop_paired ? 0x2 : 0;
		flag |= aln->secondary_aln ? 0x100 : 0; // secondary alignment
	}

	// supplementary alignment -- i.e., additional alignments that are not marked as secondary
	flag |= (i_aln > 0 && !aln->secondary_aln) ? 0x800 : 0;

	kputs(read->id, output); kputc('\t', output); // QNAME\t
	kputw((flag & 0xffff), output); kputc('\t', output); // FLAG

	if (aln->contig) { // with coordinate
		kputs(aln->contig->name, output); kputc('\t', output); // RNAME
		kputl(aln->pos, output); kputc('\t', output); // POS
		kputw(aln->mapq, output); kputc('\t', output); // MAPQ
		// BWA forces hard clipping for supplementary alignments -- i.e., additional
		// alignments that are not marked as secondary.  Those alignments are have the bit 0x800
		rapi_put_cigar(aln->n_cigar_ops, aln->cigar_ops, (i_aln > 0 && !aln->secondary_aln) ? 1 : 0, output);
	}
	else
		kputsn("*\t0\t0\t*", 7, output); // unmapped

	kputc('\t', output);

	// print the mate chr, position, and isize if applicable
	if (mate_aln->contig) {
		if (aln->contig == mate_aln->contig)
			kputc('=', output);
		else
			kputs(mate_aln->contig->name, output); // RNAME
		kputc('\t', output);
		kputl(mate_aln->pos, output); kputc('\t', output); // mate pos

		if (aln->mapped && (aln->contig == mate_aln->contig))
			kputl(rapi_get_insert_size(aln, mate_aln), output);
		else
			kputc('0', output);
	}
	else
		kputsn("*\t0\t0", 5, output);
	kputc('\t', output);

	// print SEQ and QUAL
	if (aln->secondary_aln) { // for secondary alignments, don't write SEQ and QUAL
		kputsn("*\t*", 3, output);
	}
	else {
		int i, end = read->length;
		int front_trim = 0, rear_trim = 0;
		// Trim the printed sequence for supplementary alignments (those after the
		// first in the list, so i_aln > 0, and not labeled as secondary 0x100)
		if (aln->n_cigar_ops > 0) {
			if (i_aln > 0 && (aln->cigar_ops[0].op == RAPI_CIG_S || aln->cigar_ops[0].op == RAPI_CIG_H)) {
				front_trim = aln->cigar_ops[0].len;
			}
			if (i_aln > 0 && (aln->cigar_ops[aln->n_cigar_ops - 1].op == RAPI_CIG_S || aln->cigar_ops[aln->n_cigar_ops - 1].op == RAPI_CIG_H)) {
				rear_trim = aln->cigar_ops[aln->n_cigar_ops - 1].len;
			}
		}
		int trimmed_length = read->length - front_trim - rear_trim;
		int new_size = output->l + trimmed_length + 1; // +1 for delimiter
		if (read->qual)
			new_size += trimmed_length + 1; // more room for the qual sequence
		ks_resize(output, new_size);

		if (!aln->reverse_strand) { // the forward strand
			// forward strand is simple:  front and rear trimming done to natural
			// start and end of the sequence.
			kputsn(read->seq + front_trim, trimmed_length, output);
			kputc('\t', output);
			if (read->qual) { // print qual
				for (i = front_trim; i < end - rear_trim; ++i) output->s[output->l++] = read->qual[i];
				output->s[output->l] = 0;
			}
			else kputc('*', output);
		}
		else { // the reverse strand
			// For reads on reverse strand, the CIGAR is applied backwards with respect to
			// the read->seq array.  The front_trim is applied to the end of the read and the
			// rear_trim is applied to the start.  Moreover, we have to print the reverse complement
			// of the read, so we "print" the bases in reverse order, while complementing by
			// indexing into the TGCAN char array.
			for (i = end - front_trim - 1; i >= rear_trim; --i) output->s[output->l++] = "TGCAN"[nst_nt4_table[(int)read->seq[i]]];
			kputc('\t', output);
			if (read->qual) { // print qual
				for (i = end - front_trim - 1; i >= rear_trim; --i) output->s[output->l++] = read->qual[i];
				output->s[output->l] = 0;
			}
			else kputc('*', output);
		}
	}

	// print optional tags
	if (aln->n_cigar_ops > 0) {
		kputsn("\tNM:i:", 6, output); kputw(aln->n_mismatches, output);
	}

	if (aln->score >= 0) { kputsn("\tAS:i:", 6, output); kputw(aln->score, output); }

	rapi_error_t error = RAPI_NO_ERROR;

	// write all othere tags
	for (int t = 0; t < kv_size(aln->tags) && RAPI_NO_ERROR == error; ++t) {
		kputc('\t', output);
		error = rapi_format_tag(&kv_A(aln->tags, t), output);
	}

	return error;
}

/*
 * Format SAM for a single read, using the first alignment in the
 * rapi_read->alignments list.
 *
 * \param read_num Refers to `read`. Should be 1 or 2.
 */
static rapi_error_t _rapi_format_sam_read(const rapi_read* read, const rapi_read* mate, int read_num, kstring_t* output)
{
	if (NULL == read) {
		PERROR("_rapi_format_sam_read: NULL read pointer\n");
		return RAPI_PARAM_ERROR;
	}
	rapi_error_t error = RAPI_NO_ERROR;

	if (read->n_alignments == 0) {
		error = _rapi_format_sam_aln(read, -1, mate, read_num, output);
	}
	else {
		for (int i = 0; i < read->n_alignments && !error; ++i) {
			if (i > 0) kputc('\n', output);
			error = _rapi_format_sam_aln(read, i, mate, read_num, output);
		}
	}

	return error;
}


/*
 * Format SAM for an entire fragment.
 *
 * SAM read records contain information that depends on other reads in the
 * same template (e.g., insert size, alignment coordinates for other reads,
 * flags for first/last read in template, etc.).  Generating all same for an
 * entire template within the context of a single function call makes it feasible
 * to implement this without changing the API (all the necessary info should
 * already be in here).
 *
 * However, BWA currently supports single and paired reads, so that's all we're
 * implementing in this function.
 *
 * \param reads: pointer to list of reads; reads must be ordered first to last
 * \param n_reads: number of reads in list; MUST be 1 or 2
 * \param output: str where output will be written
 */

rapi_error_t rapi_format_sam(const rapi_read** reads, int n_reads, kstring_t* output)
{
	if (n_reads <= 0 || n_reads > 2) {
		PERROR("n_reads must be 1 or 2\n");
		return RAPI_PARAM_ERROR;
	}

	rapi_error_t error = _rapi_format_sam_read(reads[0], (n_reads == 1 ? NULL : reads[1]), 1, output);
	if (n_reads == 2 && RAPI_NO_ERROR == error) {
		kputc('\n', output);
		error = _rapi_format_sam_read(reads[1], reads[0], 2, output);
	}
	return error;
}

/**
 * Format SAM for an entire fragment.
 *
 * SAM read records contain information that depends on other reads in the
 * same template (e.g., insert size, alignment coordinates for other reads,
 * flags for first/last read in template, etc.).  Generating all same for an
 * entire template within the context of a single function call makes it feasible
 * to implement this without changing the API (all the necessary info should
 * already be in here).
 *
 * However, BWA currently supports single and paired reads, so that's all we're
 * implementing in this function.
 */
rapi_error_t rapi_format_sam_b(const rapi_batch* batch, rapi_ssize_t n_frag, kstring_t* output)
{
	///// validate function arguments
	if (NULL == batch || NULL == output) {
		PERROR("NULL argument!\n");
		return RAPI_PARAM_ERROR;
	}

	if (batch->n_reads_frag > 2 || batch->n_reads_frag <= 0) {
		PERROR("Only single and paired reads are supported (got %d)\n", batch->n_reads_frag);
		return RAPI_PARAM_ERROR;
	}

	//// get read and mate

	int i_read = 0, i_mate = 1;
	const int n_reads = batch->n_reads_frag;
	const rapi_read* reads[2] = { NULL, NULL };

	reads[0] = rapi_get_read(batch, n_frag, i_read);
	if (n_reads > 1)
		reads[1] = rapi_get_read(batch, n_frag, i_mate);

	// check for errors retrieving reads
	if (NULL == reads[0] || (n_reads > 1 && NULL == reads[1])) {
		PERROR("Error fetching reads for fragment %lld: read is %s NULL; mate is %s NULL. Batch n_reads_frag: %d; n_frags: %lld.\n",
		        n_frag, (reads[0] != NULL ? "not" : ""), (reads[1] != NULL ? "not" : ""),
		        n_reads, batch->n_frags);
		return RAPI_GENERIC_ERROR;
	}

	rapi_error_t error = _rapi_format_sam_read(reads[0], (n_reads == 1 ? NULL : reads[1]), 1, output);
	if (n_reads == 2 && RAPI_NO_ERROR == error) {
		kputc('\n', output);
		error = _rapi_format_sam_read(reads[1], reads[0], 2, output);
	}
	return error;
}


rapi_error_t rapi_format_sam_hdr(const rapi_ref* ref, kstring_t* output)
{
	for (int i = 0; i < ref->n_contigs; ++i)
		ksprintf(output, "@SQ\tSN:%s\tLN:%lld\n", ref->contigs[i].name, ref->contigs[i].len);

	ksprintf(output, "@PG\tID:rapi (%s)\tPN:rapi (%s)\tVN:%s (%s)\n",
	    rapi_aligner_name(),
	    rapi_aligner_name(),
	    rapi_plugin_version(),
	    rapi_aligner_version());
	kputs("@CO File generated through the RAPI aligner interface using the specified aligner plug-in", output);
	return RAPI_NO_ERROR;
}

/**********************************/


/******** Internal structures *****/
typedef struct {
	rapi_ssize_t n_bases;
	rapi_ssize_t n_reads;
	int n_reads_per_frag;
	bseq1_t* seqs;
} bwa_batch;

static void _print_bwa_batch(FILE* out, const bwa_batch* read_batch)
{
	fprintf(out, "batch with %lld bases, %lld reads, %d reads per fragment",
	        read_batch->n_bases, read_batch->n_reads, read_batch->n_reads_per_frag);
	if (read_batch->n_reads_per_frag > 0)
		fprintf(out, " (so %lld fragments)\n", read_batch->n_reads / read_batch->n_reads_per_frag);
	else
		fprintf(out, "\n");

	fprintf(out, "=== Reads: ===\n");
	for (rapi_ssize_t r = 0; r < read_batch->n_reads; ++r) {
		const bseq1_t* bwa_read = read_batch->seqs + r;
		fprintf(out, "-=-=--=\n");
		fprintf(out, "name: %s\n", bwa_read->name);
		fprintf(out, "seq: %s\n", bwa_read->seq);
		fprintf(out, "qual %.*s\n", bwa_read->l_seq, bwa_read->qual);
	}
}

/**********************************/

/**
 * Definition of the aligner state structure.
 */
struct rapi_aligner_state {
	const library_opts* opts;
	int64_t n_reads_processed;
	// paired-end stats
	mem_pestat_t pes[4];
};


#if 1 // strdup is not defined if we compile with c99, but is when I include kvec.h.
// Is there a preprocessor DEFINE I can use to test whether a function is defined?
char* strdup(const char* str)
{
	if (NULL == str)
		return NULL;

	size_t len = strlen(str);
	char* new_str = malloc(len + 1);
	if (NULL == new_str)
		return NULL;
	return strcpy(new_str, str);
}
#endif

static rapi_error_t _library_opts_init(void) {
    _g_library_opts = calloc(1, sizeof(library_opts));
    if (!_g_library_opts) {
        PERROR("_library_opts_init Failed to allocate space for library options");
        return RAPI_MEMORY_ERROR;
    }
    return RAPI_NO_ERROR;
}

static rapi_error_t _library_opts_free(void) {
    if (_g_library_opts) {
        if (_g_library_opts->bwa_opts)
            free(_g_library_opts->bwa_opts);
        free(_g_library_opts);
        _g_library_opts = NULL;
    }
    return RAPI_NO_ERROR;
}

static rapi_error_t _set_library_opts(library_opts* lib_opts, const rapi_opts* opts) {
	lib_opts->mapq_min = opts->mapq_min;
	lib_opts->isize_min = opts->isize_min;
	lib_opts->isize_max = opts->isize_max;
	lib_opts->n_threads = opts->n_threads;
	lib_opts->bwa_opts = mem_opt_init();
	if (NULL == lib_opts->bwa_opts)
		return RAPI_MEMORY_ERROR;

	if (opts->_private) {
		// XXX: copying the mem_opt_t structure is simple at the moment,
		// but keep an eye out on future changes to it.
		const mem_opt_t* src = (mem_opt_t*)opts->_private;
		memcpy(lib_opts->bwa_opts, src, sizeof(mem_opt_t));
	}

	return RAPI_NO_ERROR;
}

static inline const library_opts* _library_opts_get(void) {
	return _g_library_opts;
}

/* Init Library */
rapi_error_t rapi_init(const rapi_opts* opts)
{
	_library_opts_free();
	rapi_error_t error = _library_opts_init();
	if (RAPI_NO_ERROR != error)
		return error;

	if (opts)
		return _set_library_opts(_g_library_opts, opts);
	else
		return RAPI_NO_ERROR;
}

rapi_error_t rapi_shutdown(void) {
	_library_opts_free();

	return RAPI_NO_ERROR;
}

/* Init Library Options */
rapi_error_t rapi_opts_init( rapi_opts * my_opts )
{
	// create a BWA opt structure
	mem_opt_t*const bwa_opt = mem_opt_init();
	if (NULL == bwa_opt)
		return RAPI_MEMORY_ERROR;

	// Default values copied from bwamem.c in 0.7.8
	bwa_opt->flag = 0;
	bwa_opt->a = 1; bwa_opt->b = 4;
	bwa_opt->o_del = bwa_opt->o_ins = 6;
	bwa_opt->e_del = bwa_opt->e_ins = 1;
	bwa_opt->w = 100;
	bwa_opt->T = 30;
	bwa_opt->zdrop = 100;
	bwa_opt->pen_unpaired = 17;
	bwa_opt->pen_clip5 = bwa_opt->pen_clip3 = 5;
	bwa_opt->min_seed_len = 19;
	bwa_opt->split_width = 10;
	bwa_opt->max_occ = 10000;
	bwa_opt->max_chain_gap = 10000;
	bwa_opt->max_ins = 10000;
	bwa_opt->mask_level = 0.50;
	bwa_opt->chain_drop_ratio = 0.50;
	bwa_opt->split_factor = 1.5;
	bwa_opt->chunk_size = 10000000;
	bwa_opt->n_threads = 1;
	bwa_opt->max_matesw = 100;
	bwa_opt->mask_level_redun = 0.95;
	bwa_opt->mapQ_coef_len = 50; bwa_opt->mapQ_coef_fac = log(bwa_opt->mapQ_coef_len);
	bwa_fill_scmat(bwa_opt->a, bwa_opt->b, bwa_opt->mat);

	my_opts->_private = bwa_opt; // hand the bwa structure onto the external one
	my_opts->ignore_unsupported = 1;
	my_opts->mapq_min     = 0;
	my_opts->isize_min    = 0;
	my_opts->isize_max    = bwa_opt->max_ins;
	my_opts->n_threads    = 1;
	kv_init(my_opts->parameters);

	return RAPI_NO_ERROR;
}

rapi_error_t rapi_opts_free( rapi_opts * my_opts )
{
	free(my_opts->_private);
	return RAPI_NO_ERROR;
}

const char* rapi_aligner_name(void)
{
	return "bwa-mem";
}

const char* rapi_aligner_version(void)
{
	return WRAPPED_BWA_VERSION;
}

const char* rapi_plugin_version(void)
{
	return RAPI_BWA_PLUGIN_VERSION;
}

/* Load Reference */
rapi_error_t rapi_ref_load( const char * reference_path, rapi_ref * ref_struct )
{
	if ( NULL == ref_struct || NULL == reference_path )
		return RAPI_PARAM_ERROR;

	const bwaidx_t*const bwa_idx = bwa_idx_load(reference_path, BWA_IDX_ALL, 1 /* use_mmap */);
	if ( NULL == bwa_idx )
		return RAPI_GENERIC_ERROR;

	// allocate memory
	ref_struct->path = strdup(reference_path);
	ref_struct->contigs = calloc( bwa_idx->bns->n_seqs, sizeof(rapi_contig) );
	if ( NULL == ref_struct->path || NULL == ref_struct->contigs )
	{
		// if either allocations we free everything and return an error
		bwa_idx_destroy((bwaidx_t*)bwa_idx);
		free(ref_struct->path);
		free(ref_struct->contigs);
		memset(ref_struct, 0, sizeof(*ref_struct));
		return RAPI_MEMORY_ERROR;
	}

	/* Fill in Contig Information */
	ref_struct->n_contigs = bwa_idx->bns->n_seqs; /* Handle contains bnt_seq_t * bns holding contig information */
	for ( int i = 0; i < ref_struct->n_contigs; ++i )
	{
		rapi_contig* c = &ref_struct->contigs[i];
		c->len = bwa_idx->bns->anns[i].len;
		c->name = bwa_idx->bns->anns[i].name; // points to BWA string
		c->assembly_identifier = NULL;
		c->species = NULL;
		c->uri = NULL;
		c->md5 = NULL;
	}
	ref_struct->_private = (bwaidx_t*)bwa_idx;

	return RAPI_NO_ERROR;
}

/* Free Reference */
rapi_error_t rapi_ref_free( rapi_ref * ref )
{
	// free bwa's part
	bwa_idx_destroy(ref->_private);

	// then free the rest of the structure
	free(ref->path);

	if (ref->contigs) {
		for ( int i = 0; i < ref->n_contigs; ++i )
		{
			rapi_contig* c = &ref->contigs[i];
			// *Don't* free name since it points to BWA's string
			free(c->assembly_identifier);
			free(c->species);
			free(c->uri);
			free(c->md5);
		}
		free (ref->contigs);
	}
	memset(ref, 0, sizeof(*ref));
	return RAPI_NO_ERROR;
}

void _free_bwa_batch_contents(bwa_batch* batch)
{
	for (int i = 0; i < batch->n_reads; ++i) {
		// *Don't* free the read id since it points to the original structure's string
		free(batch->seqs[i].seq);
		free(batch->seqs[i].qual);
	}

	free(batch->seqs);
	memset(batch, 0, sizeof(bwa_batch));
}


static rapi_error_t _batch_to_bwa_seq(const rapi_batch* batch, const rapi_opts* opts, int start_fragment, int end_fragment, bwa_batch* bwa_seqs)
{
	if (start_fragment < 0 && end_fragment < 0) {
		start_fragment = 0;
		end_fragment = batch->n_frags;
	}

	if (end_fragment > batch->n_frags || start_fragment > end_fragment) {
		PERROR("start or end fragmet is out of bounds. Got start %d and end %d but we have %lld fragments\n",
		        start_fragment, end_fragment, batch->n_frags);
		return RAPI_PARAM_ERROR;
	}

	bwa_seqs->n_bases = 0;
	bwa_seqs->n_reads = 0;
	bwa_seqs->n_reads_per_frag = batch->n_reads_frag;

	int n_frags = end_fragment - start_fragment;

	bwa_seqs->seqs = calloc(n_frags * batch->n_reads_frag, sizeof(bseq1_t));
	if (NULL == bwa_seqs->seqs) {
		PERROR("Allocation failed!\n");
		return RAPI_MEMORY_ERROR;
	}

	for (int f = start_fragment; f < end_fragment; ++f)
	{
		for (int r = 0; r < batch->n_reads_frag; ++r)
		{
			const rapi_read*const rapi_read = rapi_get_read(batch, f, r);
			bseq1_t* bwa_read = bwa_seqs->seqs + bwa_seqs->n_reads;

			// -- In bseq1_t, all strings are null-terminated.
			// We duplicated the seq and qual since BWA modifies them.
			bwa_read->seq = strdup(rapi_read->seq);
			bwa_read->qual = (NULL == rapi_read->qual) ? NULL : strdup(rapi_read->qual);
			if (!bwa_read->seq || (rapi_read->qual && !bwa_read->qual))
				goto failed_allocation;

			bwa_read->name = rapi_read->id;
			bwa_read->l_seq = rapi_read->length;
			// Since we use calloc to allocate this structures there's no need to set
			// these members to NULL.
			//bwa_read->comment = NULL;
			//bwa_read->sam = NULL;

			// As we build the objects n_reads keeps the actual number
			// constructed thus far and thus can be used to free the allocated
			// structres in case of error.
			bwa_seqs->n_reads += 1;
			bwa_seqs->n_bases += rapi_read->length;
		}
	}
	return RAPI_NO_ERROR;

failed_allocation:
	PERROR("Failed to allocate while constructing sequences! Freeing and returning\n");
	_free_bwa_batch_contents(bwa_seqs);
	return RAPI_MEMORY_ERROR;
}

/*
 * Many of the default option values need to be adjusted if the matching score
 * (opt->a) is changed.  This function (from the BWA code) does that.
 *
 * \param opt The structure containing actual values and that will be "fixed"
 * \param override set the members to 1 if the value has been overridden and thus should be kept;
 *                 else the corresponding value in opts is assumed to be at default and will be adjusted.
 *
 * So, if you change 'a' set the 'override' accordingly and call this function.  E.g.,
 *
 * <pre>
 *   rapi_opts* opt = rapi_init_opts();
 *   mem_opt_t* bwa_opts = (mem_opt_t*) opt->_private;
 *   bwa_opts->a = 2;
 *   bwa_opts->b = 5;
 *   mem_opt_t override;
 *   override.a = override.b = 1;
 *   adjust_bwa_opts(bwa_opts, &override);
 * </pre>
 */
void adjust_bwa_opts(mem_opt_t* opt, const mem_opt_t* override)
{
	if (override->a == 1) { // matching score is changed
		if (override->b != 1) opt->b *= opt->a;
		if (override->T != 1) opt->T *= opt->a;
		if (override->o_del != 1) opt->o_del *= opt->a;
		if (override->e_del != 1) opt->e_del *= opt->a;
		if (override->o_ins != 1) opt->o_ins *= opt->a;
		if (override->e_ins != 1) opt->e_ins *= opt->a;
		if (override->zdrop != 1) opt->zdrop *= opt->a;
		if (override->pen_clip5 != 1) opt->pen_clip5 *= opt->a;
		if (override->pen_clip3 != 1) opt->pen_clip3 *= opt->a;
		if (override->pen_unpaired != 1) opt->pen_unpaired *= opt->a;
		bwa_fill_scmat(opt->a, opt->b, opt->mat);
	}
}

static int _convert_opts(const library_opts* opts, mem_opt_t* bwa_opts)
{
	bwa_opts->T = opts->mapq_min;
	bwa_opts->max_ins = opts->isize_max;
	bwa_opts->n_threads = opts->n_threads;

	// TODO: other options provided through 'parameters' field
	return RAPI_NO_ERROR;
}

rapi_error_t rapi_aligner_state_init(struct rapi_aligner_state** ret_state, const rapi_opts* opts)
{
	rapi_error_t error;
	const library_opts* lib_opts;
	library_opts* tmp = NULL;

	if (opts) {
		tmp = calloc(1, sizeof(library_opts));
		if (!tmp) return RAPI_MEMORY_ERROR;

		error = _set_library_opts(tmp, opts);

		if (error != RAPI_NO_ERROR) {
			free(tmp);
			return error;
		}
		lib_opts = tmp;
	}
	else {
		lib_opts = _library_opts_get();
	}

	// allocate and zero the structure
	rapi_aligner_state* state = *ret_state = calloc(1, sizeof(rapi_aligner_state));
	if (NULL == state) {
		if (_library_opts_get() != lib_opts) // if it's not the library-wide option structure
			free(tmp); // free through tmp avoid the warning about dropping 'const'
		return RAPI_MEMORY_ERROR;
	}

	state->opts = lib_opts;

	return RAPI_NO_ERROR;
}

rapi_error_t rapi_aligner_state_free(rapi_aligner_state* state)
{
	if (state->opts != _library_opts_get()) {
		free((library_opts*)state->opts);
		state->opts = NULL;
	}
	free(state);
	return RAPI_NO_ERROR;
}

void rapi_put_cigar(int n_ops, const rapi_cigar* ops, int force_hard_clip, kstring_t* output)
{
	if (n_ops > 0) {
		for (int i = 0; i < n_ops; ++i) {
			int c = ops[i].op;
			if (c == RAPI_CIG_S || c == RAPI_CIG_H) c = force_hard_clip ? RAPI_CIG_H : RAPI_CIG_S;
			kputw(ops[i].len, output);
			kputc(rapi_cigops_char[c], output);
		}
	}
	else
		kputc('*', output);
}

long rapi_get_insert_size(const rapi_alignment* read, const rapi_alignment* mate)
{
	long isize = 0;

	if (read->mapped && mate->mapped && (read->contig == mate->contig))
	{
		if (mate->n_cigar_ops == 0 || read->n_cigar_ops == 0)
			err_fatal(__func__, "No cigar ops for mapped reads! aln->n_cigar_ops: %d; mate_aln->n_cigar_ops: %d\n", read->n_cigar_ops, mate->n_cigar_ops);

		int64_t p0 = read->pos + (read->reverse_strand ? rapi_get_rlen(read->n_cigar_ops, read->cigar_ops) - 1 : 0);
		int64_t p1 = mate->pos + (mate->reverse_strand ? rapi_get_rlen(mate->n_cigar_ops, mate->cigar_ops) - 1 : 0);
		isize = -(p0 - p1 + (p0 > p1? 1 : p0 < p1? -1 : 0));
	}
	return isize;
}

int rapi_get_rlen(int n_cigar, const rapi_cigar* cigar_ops)
{
	int len = 0;
	for (int k = 0; k < n_cigar; ++k) {
		int op = cigar_ops[k].op;
		if (op == RAPI_CIG_M || op == RAPI_CIG_D)
			len += cigar_ops[k].len;
	}
	return len;
}


/********** modified BWA code *****************/

/* Copied directly from bwamem_pair */
inline int mem_infer_dir(int64_t l_pac, int64_t b1, int64_t b2, int64_t *dist)
{
	int64_t p2;
	int r1 = (b1 >= l_pac), r2 = (b2 >= l_pac);
	p2 = r1 == r2? b2 : (l_pac<<1) - 1 - b2; // p2 is the coordinate of read 2 on the read 1 strand
	*dist = p2 > b1? p2 - b1 : b1 - p2;
	return (r1 == r2? 0 : 1) ^ (p2 > b1? 0 : 3);
}

// IMPORTANT: must run mem_sort_and_dedup() before calling the mem_mark_primary_se function (but it's called by mem_align1_core)

/* based on mem_aln2sam */
static int _bwa_aln_to_rapi_aln(const rapi_ref* rapi_ref, rapi_read* our_read, int is_paired,
		const bseq1_t *s,
		const mem_aln_t *const bwa_aln_list, int list_length)
{
	if (list_length < 0)
		return RAPI_PARAM_ERROR;

	rapi_tag* pTag; // temporary pointer to form tags

	our_read->alignments = calloc(list_length, sizeof(rapi_alignment));
	if (NULL == our_read->alignments)
		return RAPI_MEMORY_ERROR;
	our_read->n_alignments = list_length;

	for (int which = 0; which < list_length; ++which)
	{
		const mem_aln_t* bwa_aln = &bwa_aln_list[which];
		rapi_alignment* our_aln = &our_read->alignments[which];

		if (bwa_aln->rid >= rapi_ref->n_contigs) { // huh?? Out of bounds
			PERROR("read reference id value %d is out of bounds (n_contigs: %d)\n", bwa_aln->rid, rapi_ref->n_contigs);
			free(our_read->alignments);
			our_read->alignments = NULL; our_read->n_alignments = 0;
			return RAPI_GENERIC_ERROR;
		}

		// set flags
		our_aln->paired = is_paired != 0;
		our_aln->prop_paired = (bwa_aln->flag & 0x2) != 0; // 0x2 is the SAM proper pair flag
		our_aln->score = bwa_aln->score;
		our_aln->mapq = bwa_aln->mapq;
		// In BWA's code (e.g., mem_aln2sam) when the 0x10000 bit is set the alignment
		// is printed as a secondary alignment (i.e., the 0x100 bit is set in the flag).
		our_aln->secondary_aln = ((bwa_aln->flag & 0x100) | (bwa_aln->flag & 0x10000)) != 0;

		our_aln->mapped = bwa_aln->rid >= 0;
		if (bwa_aln->rid >= 0) { // with coordinate
			our_aln->reverse_strand = bwa_aln->is_rev != 0;
			our_aln->contig = &rapi_ref->contigs[bwa_aln->rid];
			our_aln->pos = bwa_aln->pos + 1;
			our_aln->n_mismatches = bwa_aln->NM;
			if (bwa_aln->n_cigar) { // aligned
				our_aln->cigar_ops = malloc(bwa_aln->n_cigar * sizeof(our_aln->cigar_ops[0]));
				if (NULL == our_aln->cigar_ops)
					err_fatal(__func__, "Failed to allocate cigar space");
				our_aln->n_cigar_ops = bwa_aln->n_cigar;
				for (int i = 0; i < bwa_aln->n_cigar; ++i) {
					our_aln->cigar_ops[i].op = bwa_aln->cigar[i] & 0xf;
					our_aln->cigar_ops[i].len = bwa_aln->cigar[i] >> 4;
				}

				// BWA stores the MD string right after the cigar array.
				const char* md = (char*)(bwa_aln->cigar + bwa_aln->n_cigar);
				pTag = kv_pushp(rapi_tag, our_aln->tags);
				rapi_tag_set_key(pTag, "MD");
				rapi_tag_set_text(pTag, md);
			}
		}

		if (bwa_aln->sub >= 0) {
			pTag = kv_pushp(rapi_tag, our_aln->tags);
			rapi_tag_set_key(pTag, "XS");
			rapi_tag_set_long(pTag, bwa_aln->sub);
		}
	}

	// generate SA tags for all alignments in the list
	kstring_t tmp_sa = { 0, 0, NULL };
	for (int i_aln = 0; i_aln < our_read->n_alignments; ++i_aln) {
		tmp_sa.l = 0; // reposition the string cursor to the beginning
		rapi_alignment*const aln = our_read->alignments + i_aln;

		if (!(aln->secondary_aln)) { // not multi-hit --
			// XXX: actually, secondary_aln is set if BWA set either 0x100 or 0x10000
			// are set in the alignments' flag.  On the other hand, BWA's original code only
			// checks the 0x100 bit :-o

			// Find the first alignment in the list, which is not the main alignment that
			// we're considering (i.e., i_aln), and is not a secondary alignment
			int other_primary = 0;
			for (other_primary = 0; other_primary < our_read->n_alignments; ++other_primary) {
				if (other_primary != i_aln && !our_read->alignments[other_primary].secondary_aln)
					break;
			}
			if (other_primary < our_read->n_alignments) { // there are other primary hits; output them
				for (; other_primary < our_read->n_alignments; ++other_primary) {
					const rapi_alignment*const sa = our_read->alignments + other_primary;

					// proceed if: 1) different from the current; 2) not shadowed multi hit
					if (other_primary == i_aln || sa->secondary_aln) continue;

					kputs(sa->contig->name, &tmp_sa); kputc(',', &tmp_sa);
					kputl(sa->pos, &tmp_sa); kputc(',', &tmp_sa); // XXX: BWA has sa->pos + 1
					kputc(sa->reverse_strand ? '-' : '+', &tmp_sa); kputc(',', &tmp_sa);
					rapi_put_cigar(sa->n_cigar_ops, sa->cigar_ops, 0, &tmp_sa);
					kputc(',', &tmp_sa); kputw(sa->mapq, &tmp_sa);
					kputc(',', &tmp_sa); kputw(sa->n_mismatches, &tmp_sa);
					kputc(';', &tmp_sa);
				}

				// now set the tag
				rapi_tag* pTag = kv_pushp(rapi_tag, aln->tags);
				rapi_tag_set_key(pTag, "SA");
				rapi_tag_set_text(pTag, tmp_sa.s);
			}
		}
	}
	free(tmp_sa.s);

	return RAPI_NO_ERROR;
}

/*
 * Based on mem_reg2sam_se.
 * We took out the call to mem_aln2sam and instead write the result to
 * the corresponding rapi_read structure.
 */
static int _bwa_reg2_rapi_aln(const mem_opt_t *opt, const rapi_ref* rapi_ref, rapi_read* our_read, int is_paired, bseq1_t *seq, mem_alnreg_v *a, int extra_flag)
{
	rapi_error_t error = RAPI_NO_ERROR;
	const bntseq_t *const bns = ((bwaidx_t*)rapi_ref->_private)->bns;
	const uint8_t *const pac = ((bwaidx_t*)rapi_ref->_private)->pac;

	kvec_t(mem_aln_t) aa;
	int k;

	kv_init(aa);
	for (k = 0; k < a->n; ++k) {
		mem_alnreg_t *p = &a->a[k];
		mem_aln_t *q;
		if (p->score < opt->T) continue;
		if (p->secondary >= 0 && !(opt->flag&MEM_F_ALL)) continue;
		if (p->secondary >= 0 && p->score < a->a[p->secondary].score * .5) continue;
		q = kv_pushp(mem_aln_t, aa);
		*q = mem_reg2aln(opt, bns, pac, seq->l_seq, seq->seq, p);
		q->flag |= (is_paired ? 0x1 : 0);
		q->flag |= extra_flag; // flag secondary
		if (p->secondary >= 0) q->sub = -1; // don't output sub-optimal score
		if (k && p->secondary < 0) // if supplementary
			q->flag |= (opt->flag&MEM_F_NO_MULTI)? 0x10000 : 0x800;
		if (k && q->mapq > aa.a[0].mapq) q->mapq = aa.a[0].mapq;
	}
	if (aa.n == 0) { // no alignments good enough; then write an unaligned record
		mem_aln_t t;
		t = mem_reg2aln(opt, bns, pac, seq->l_seq, seq->seq, 0);
		t.flag |= extra_flag;
		// RAPI
		error = _bwa_aln_to_rapi_aln(rapi_ref, our_read, is_paired, seq, &t, 1);
	}
	else {
		error = _bwa_aln_to_rapi_aln(rapi_ref, our_read, is_paired, seq, /* list of aln */ aa.a, aa.n);
	}

	if (aa.n > 0)
	{
		for (int k = 0; k < aa.n; ++k)
			free(aa.a[k].cigar);
		free(aa.a);
	}
	return error;
}


#if 1

#define raw_mapq(diff, a) ((int)(6.02 * (diff) / (a) + .499))
/*
 * Mostly taken from mem_sam_pe in bwamem_pair.c
 *
 * \return I think this function returns the number pairs aligned by SW
 */
int _bwa_mem_pe(const mem_opt_t *opt, const rapi_ref* rapi_ref, const mem_pestat_t pes[4], uint64_t id, bseq1_t s[2], mem_alnreg_v a[2], rapi_read out[2])
{
	const bntseq_t *const bns = ((bwaidx_t*)rapi_ref->_private)->bns;
	const uint8_t *const pac = ((bwaidx_t*)rapi_ref->_private)->pac;

	int n = 0, i, j, z[2], o, subo, n_sub, extra_flag = 1;
	kstring_t str;
	mem_aln_t h[2];

	str.l = str.m = 0; str.s = 0;
	if (!(opt->flag & MEM_F_NO_RESCUE)) { // then perform SW for the best alignment
		mem_alnreg_v b[2];
		kv_init(b[0]); kv_init(b[1]);
		for (i = 0; i < 2; ++i)
			for (j = 0; j < a[i].n; ++j)
				if (a[i].a[j].score >= a[i].a[0].score  - opt->pen_unpaired)
					kv_push(mem_alnreg_t, b[i], a[i].a[j]);
		for (i = 0; i < 2; ++i)
			for (j = 0; j < b[i].n && j < opt->max_matesw; ++j)
				n += mem_matesw(opt, bns->l_pac, pac, pes, &b[i].a[j], s[!i].l_seq, (uint8_t*)s[!i].seq, &a[!i]);
		free(b[0].a); free(b[1].a);
	}
	mem_mark_primary_se(opt, a[0].n, a[0].a, id<<1|0);
	mem_mark_primary_se(opt, a[1].n, a[1].a, id<<1|1);
	if (opt->flag&MEM_F_NOPAIRING) goto no_pairing;
	// pairing single-end hits
	if (a[0].n && a[1].n && (o = mem_pair(opt, bns->l_pac, pac, pes, s, a, id, &subo, &n_sub, z)) > 0) {
		int is_multi[2], q_pe, score_un, q_se[2];
		// check if an end has multiple hits even after mate-SW
		for (i = 0; i < 2; ++i) {
			for (j = 1; j < a[i].n; ++j)
				if (a[i].a[j].secondary < 0 && a[i].a[j].score >= opt->T) break;
			is_multi[i] = j < a[i].n? 1 : 0;
		}
		if (is_multi[0] || is_multi[1]) goto no_pairing; // TODO: in rare cases, the true hit may be long but with low score
		// compute mapQ for the best SE hit
		score_un = a[0].a[0].score + a[1].a[0].score - opt->pen_unpaired;
		//q_pe = o && subo < o? (int)(MEM_MAPQ_COEF * (1. - (double)subo / o) * log(a[0].a[z[0]].seedcov + a[1].a[z[1]].seedcov) + .499) : 0;
		subo = subo > score_un? subo : score_un;
		q_pe = raw_mapq(o - subo, opt->a);
		if (n_sub > 0) q_pe -= (int)(4.343 * log(n_sub+1) + .499);
		if (q_pe < 0) q_pe = 0;
		if (q_pe > 60) q_pe = 60;
		// the following assumes no split hits
		if (o > score_un) { // paired alignment is preferred
			mem_alnreg_t *c[2];
			c[0] = &a[0].a[z[0]]; c[1] = &a[1].a[z[1]];
			for (i = 0; i < 2; ++i) {
				if (c[i]->secondary >= 0)
					c[i]->sub = a[i].a[c[i]->secondary].score, c[i]->secondary = -2;
				q_se[i] = mem_approx_mapq_se(opt, c[i]);
			}
			q_se[0] = q_se[0] > q_pe? q_se[0] : q_pe < q_se[0] + 40? q_pe : q_se[0] + 40;
			q_se[1] = q_se[1] > q_pe? q_se[1] : q_pe < q_se[1] + 40? q_pe : q_se[1] + 40;
			extra_flag |= 2;
			// cap at the tandem repeat score
			q_se[0] = q_se[0] < raw_mapq(c[0]->score - c[0]->csub, opt->a)? q_se[0] : raw_mapq(c[0]->score - c[0]->csub, opt->a);
			q_se[1] = q_se[1] < raw_mapq(c[1]->score - c[1]->csub, opt->a)? q_se[1] : raw_mapq(c[1]->score - c[1]->csub, opt->a);
		} else { // the unpaired alignment is preferred
			z[0] = z[1] = 0;
			q_se[0] = mem_approx_mapq_se(opt, &a[0].a[0]);
			q_se[1] = mem_approx_mapq_se(opt, &a[1].a[0]);
		}

		// write SAM
		h[0] = mem_reg2aln(opt, bns, pac, s[0].l_seq, s[0].seq, &a[0].a[z[0]]); h[0].mapq = q_se[0]; h[0].flag |= 0x40 | extra_flag;
		h[1] = mem_reg2aln(opt, bns, pac, s[1].l_seq, s[1].seq, &a[1].a[z[1]]); h[1].mapq = q_se[1]; h[1].flag |= 0x80 | extra_flag;
		// RAPI: instead of writing sam, convert mem_aln_t into our alignments
		// XXX: I'm not so sure about the alignment I'm passing in.  Review
		int error1 = _bwa_aln_to_rapi_aln(rapi_ref, &out[0], 1, &s[0], &h[0], 1);
		int error2 = _bwa_aln_to_rapi_aln(rapi_ref, &out[1], 1, &s[1], &h[1], 1);
		if (error1 || error2) {
			err_fatal(__func__, "error %d while converting BWA mem_aln_t for read %d into rapi alignments\n", (error1 ? 1 : 2), (error1 ? error1 : error2));
			abort();
		}

		if (strcmp(s[0].name, s[1].name) != 0) err_fatal(__func__, "paired reads have different names: \"%s\", \"%s\"\n", s[0].name, s[1].name);
		free(h[0].cigar); free(h[1].cigar);

	} else goto no_pairing;
	return n;

no_pairing:
	for (i = 0; i < 2; ++i) {
		if (a[i].n && a[i].a[0].score >= opt->T)
			h[i] = mem_reg2aln(opt, bns, pac, s[i].l_seq, s[i].seq, &a[i].a[0]);
		else h[i] = mem_reg2aln(opt, bns, pac, s[i].l_seq, s[i].seq, 0);
	}
	if (!(opt->flag & MEM_F_NOPAIRING) && h[0].rid == h[1].rid && h[0].rid >= 0) { // if the top hits from the two ends constitute a proper pair, flag it.
		int64_t dist;
		int d;
		d = mem_infer_dir(bns->l_pac, a[0].a[0].rb, a[1].a[0].rb, &dist);
		if (!pes[d].failed && dist >= pes[d].low && dist <= pes[d].high) extra_flag |= 2;
	}

	// We need to pass the extra flag bits to _bwa_reg2_rapi_aln because it needs to set them
	// on any secondary alignments.
	int error1 = _bwa_reg2_rapi_aln(opt, rapi_ref, &out[0], 1, &s[0], &a[0], 0x41|extra_flag);
	int error2 = _bwa_reg2_rapi_aln(opt, rapi_ref, &out[1], 1, &s[1], &a[1], 0x81|extra_flag);
	if (error1 || error2) {
		err_fatal(__func__, "error %d while converting *with no pairing* BWA mem_aln_t for read %d into rapi alignments\n", (error1 ? 1 : 2), (error1 ? error1 : error2));
		abort();
	}

	if (strcmp(s[0].name, s[1].name) != 0) err_fatal(__func__, "paired reads have different names: \"%s\", \"%s\"\n", s[0].name, s[1].name);
	free(h[0].cigar); free(h[1].cigar);
	return n;
}

typedef struct {
	const mem_opt_t *opt;
	const rapi_ref* rapi_ref;
	const bwa_batch* read_batch;
	rapi_read* rapi_reads; // need to pass these along because the code to convert BWA alignments into rapi is nested pretty deep
	mem_pestat_t *pes;
	mem_alnreg_v *regs;
	int64_t n_processed;
} bwa_worker_t;

/*
 * This function is the same as worker1 from bwamem.c
 */
static void bwa_worker_1(void *data, int i, int tid)
{
	bwa_worker_t *w = (bwa_worker_t*)data;

	const bwaidx_t* const bwaidx = (bwaidx_t*)(w->rapi_ref->_private);
	const bwt_t*    const bwt    = bwaidx->bwt;
	const bntseq_t* const bns    = bwaidx->bns;
	const uint8_t*  const pac    = bwaidx->pac;

	//PDEBUG("bwa_worker_1: MEM_F_PE is %sset\n", ((w->opt->flag & MEM_F_PE) == 0 ? "not " : " "));
	if (w->opt->flag & MEM_F_PE) {
		int read = 2*i;
		int mate = 2*i + 1;
		w->regs[read] = mem_align1_core(w->opt, bwt, bns, pac, w->read_batch->seqs[read].l_seq, w->read_batch->seqs[read].seq);
		w->regs[mate] = mem_align1_core(w->opt, bwt, bns, pac, w->read_batch->seqs[mate].l_seq, w->read_batch->seqs[mate].seq);
	} else {
		w->regs[i] = mem_align1_core(w->opt, bwt, bns, pac, w->read_batch->seqs[i].l_seq, w->read_batch->seqs[i].seq);
	}
}

/* based on worker2 from bwamem.c */
static void bwa_worker_2(void *data, int i, int tid)
{
	bwa_worker_t *w = (bwa_worker_t*)data;
	//PDEBUG("bwa_worker_2 with i %d\n", i);
	rapi_error_t error = RAPI_NO_ERROR;

	if ((w->opt->flag & MEM_F_PE)) {
		// paired end
		// This function does not return an error code but aborts if things go wrong.
		// Unfortunately this strategy is nested deep in the BWA code.
		//mem_sam_pe(w->opt, w->bns, w->pac, w->pes, (w->n_processed>>1) + i, &w->seqs[i<<1], &w->regs[i<<1]);
		_bwa_mem_pe(w->opt, w->rapi_ref, w->pes, w->n_processed / 2 + i,
		            &(w->read_batch->seqs[2 * i]), &w->regs[2 * i], &(w->rapi_reads[2 * i]));
		free(w->regs[2 * i].a); kv_init(w->regs[2 * i]);
		free(w->regs[2 * i + 1].a); kv_init(w->regs[2 * i + 1]);
	}
	else {
		// single end
		PERROR("Single end alignments aren't implemented in rapi_bwa yet!");
		error = RAPI_OP_NOT_SUPPORTED_ERROR;
		mem_mark_primary_se(w->opt, w->regs[i].n, w->regs[i].a, w->n_processed + i);
		//mem_reg2sam_se(w->opt, w->bns, w->pac, &w->seqs[i], &w->regs[i], 0, 0);
		//error = _bwa_reg2_rapi_aln(w->opt, w->rapi_ref, &(w->read_batch->seqs[i]), /* unpaired */ 0, &w->regs[i], &(w->rapi_reads[i]), 0, 0);
		free(w->regs[i].a); kv_init(w->regs[i]);
	}

	if (error != RAPI_NO_ERROR) {
		err_fatal(__func__, "%s (%d) while running %s end alignments\n",
		rapi_error_name(error), error, ((w->opt->flag & MEM_F_PE) ? "pair" : "single"));
	}
}

#endif
/********** end modified BWA code *****************/

rapi_error_t rapi_align_reads( const rapi_ref* ref, rapi_batch* batch,
        rapi_ssize_t start_fragment, rapi_ssize_t end_fragment, rapi_aligner_state* state )
{
	rapi_error_t error = RAPI_NO_ERROR;

	if (batch->n_reads_frag > 2)
		return RAPI_OP_NOT_SUPPORTED_ERROR;

	if (batch->n_reads_frag <= 0)
		return RAPI_PARAM_ERROR;

	// "extract" BWA-specific structures
	mem_opt_t*const bwa_opt = (mem_opt_t*) state->opts->bwa_opts;

	if (batch->n_reads_frag == 2) // paired-end
		bwa_opt->flag |= MEM_F_PE;

	if ((error = _convert_opts(state->opts, bwa_opt)))
		return error;

	// traslate our read structure into BWA reads
	bwa_batch bwa_seqs;
	if ((error = _batch_to_bwa_seq(batch, state->opts, start_fragment, end_fragment, &bwa_seqs)))
		return error;
	fprintf(stderr, "Converted reads to BWA structures.\n");

	fprintf(stderr, "Going to process.\n");
	mem_alnreg_v *regs = malloc(bwa_seqs.n_reads * sizeof(mem_alnreg_v));
	if (NULL == regs) {
		error = RAPI_MEMORY_ERROR;
		goto clean_up;
	}

	extern void kt_for(int n_threads, void (*func)(void*,int,int), void *data, int n);
	bwa_worker_t w;
	w.opt = bwa_opt;
	w.read_batch = &bwa_seqs;
	w.regs = regs;
	w.pes = state->pes;
	w.n_processed = state->n_reads_processed;
	w.rapi_ref = ref;
	w.rapi_reads = batch->reads;

	fprintf(stderr, "Calling bwa_worker_1. ");
	rapi_print_bwa_flag_string(stderr, bwa_opt->flag);

	int n_fragments = (bwa_opt->flag & MEM_F_PE) ? bwa_seqs.n_reads / 2 : bwa_seqs.n_reads;
	kt_for(bwa_opt->n_threads, bwa_worker_1, &w, n_fragments); // find mapping positions

	if (bwa_opt->flag & MEM_F_PE) { // infer insert sizes if not provided
		// TODO: support manually setting insert size dist parameters
		// if (pes0) memcpy(pes, pes0, 4 * sizeof(mem_pestat_t)); // if pes0 != NULL, set the insert-size distribution as pes0
		mem_pestat(bwa_opt, ((bwaidx_t*)ref->_private)->bns->l_pac, bwa_seqs.n_reads, regs, w.pes); // infer the insert size distribution from data
	}
	kt_for(bwa_opt->n_threads, bwa_worker_2, &w, n_fragments); // generate alignment

	// run the alignment
	state->n_reads_processed += bwa_seqs.n_reads;
	fprintf(stderr, "processed %" PRId64 " reads\n", state->n_reads_processed);

clean_up:
	free(regs);
	_free_bwa_batch_contents(&bwa_seqs);

	return error;
}

/******* Generic functions ***********
 * These are generally usable. We should put them in a generic rapi.c file.
 *************************************/

/* Allocate reads */
rapi_error_t rapi_reads_alloc( rapi_batch * batch, int n_reads_fragment, int n_fragments )
{
	if (n_fragments < 0 || n_reads_fragment < 0)
		return RAPI_PARAM_ERROR;

	batch->reads = calloc( n_reads_fragment * n_fragments, sizeof(rapi_read) );
	if (NULL == batch->reads)
		return RAPI_MEMORY_ERROR;
	batch->n_frags = n_fragments;
	batch->n_reads_frag = n_reads_fragment;
	return RAPI_NO_ERROR;
}

rapi_error_t rapi_reads_reserve(rapi_batch* batch, rapi_ssize_t n_fragments)
{
	if (n_fragments < 0)
		return RAPI_PARAM_ERROR;
	if (n_fragments == 0)
		return RAPI_NO_ERROR;

	if (n_fragments > batch->n_frags)
	{
		// Current space insufficient.  Need to reallocate.
		rapi_ssize_t old_n_reads = batch->n_frags * batch->n_reads_frag;
		rapi_ssize_t new_n_reads = n_fragments * batch->n_reads_frag;
		rapi_read* space = realloc(batch->reads, new_n_reads * sizeof(batch->reads[0]));
		if (space == NULL)
			return RAPI_MEMORY_ERROR;
		else {
			// set new space to 0
			memset(space + old_n_reads, 0, (new_n_reads - old_n_reads) * sizeof(batch->reads[0]));
			batch->n_frags = n_fragments;
			batch->reads = space;
		}
	}
	return RAPI_NO_ERROR;
}

static void _rapi_free_read_structures(rapi_batch* batch)
{
	for (rapi_ssize_t f = 0; f < batch->n_frags; ++f) {
		for (int r = 0; r < batch->n_reads_frag; ++r) {
			rapi_read* read = rapi_get_read(batch, f, r);
			// the reads use a single chunk of memory for id, seq and quality
			free(read->id);
			for (int a = 0; a < read->n_alignments; ++a) {
				for (int t = 0; t < read->alignments[a].tags.n; ++t)
					rapi_tag_clear(&read->alignments[a].tags.a[t]);
				kv_destroy(read->alignments[a].tags);
				free(read->alignments[a].cigar_ops);
			}
			free(read->alignments);
			read->n_alignments = 0;
			// *Don't* free the contig name.  It belongs to the contig structure.
		}
	}
}

rapi_error_t rapi_reads_clear(rapi_batch* batch)
{
	_rapi_free_read_structures(batch);
	memset(batch->reads, 0,  batch->n_reads_frag * batch->n_frags * sizeof(rapi_read));

	return RAPI_NO_ERROR;
}

rapi_error_t rapi_reads_free(rapi_batch* batch )
{
	_rapi_free_read_structures(batch);
	free(batch->reads);
	memset(batch, 0, sizeof(*batch));

	return RAPI_NO_ERROR;
}

rapi_error_t rapi_set_read(rapi_batch* batch,
	        rapi_ssize_t n_frag, int n_read,
	        const char* name, const char* seq, const char* qual,
	        int q_offset) {
	rapi_error_t error_code = RAPI_NO_ERROR;

	if (n_frag < 0 || n_frag >= batch->n_frags ||
	    n_read < 0 || n_read >= batch->n_reads_frag)
		return RAPI_PARAM_ERROR;

	rapi_read* read = rapi_get_read(batch, n_frag, n_read);
	const int name_len = strlen(name);

	const int seq_len = strlen(seq);
	if (seq_len == 0) {
		PERROR("Got sequence of length 0\n");
		return RAPI_PARAM_ERROR;
	}

	read->length = seq_len;

	// simplify allocation and error checking by allocating a single buffer
	int buf_size = name_len + 1 + seq_len + 1;
	if (qual)
		buf_size += seq_len + 1;

	read->id = malloc(buf_size);
	if (NULL == read->id) { // failed allocation
		PERROR("Unable to allocate memory for sequence\n");
		return RAPI_MEMORY_ERROR;
	}

	// copy name
	strcpy(read->id, name);

	// sequence, placed right after the name
	read->seq = read->id + name_len + 1;
	strcpy(read->seq, seq);

	// the quality, if we have it, may need to be recoded
	if (NULL == qual)
		read->qual = NULL;
	else {
		read->qual = read->seq + seq_len + 1;
		for (int i = 0; i < seq_len; ++i) {
			read->qual[i] = (int)qual[i] - q_offset + 33; // 33 is the Sanger offset.  BWA expects it this way.
			if (read->qual[i] < 33 || read->qual[i] > 126)
			{ // Sanger base qualities have an allowed range of [0,93], and 93+33=126
				PERROR("Invalid base quality score %d\n", read->qual[i]);
				error_code = RAPI_PARAM_ERROR;
				goto error;
			}
		}
		read->qual[seq_len] = '\0';
	}

	// trim from the name /[12]$
	int t = name_len;
	if (t > 2 && read->id[t-2] == '/' && (read->id[t-1] == '1' || read->id[t-1] == '2'))
		read->id[t-2] = '\0';

	return RAPI_NO_ERROR;

error:
	// In case of error, free any allocated memory and return the error
	free(read->id);
	read->id = NULL;
	return error_code;
}

