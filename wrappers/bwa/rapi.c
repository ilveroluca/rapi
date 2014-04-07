

#include "../../aligner.h"

#include <stdio.h>
#include <stdlib.h>

/* SAM flags */
#define SAM_FPD   1 // paired
#define SAM_FPP   2 // properly paired
#define SAM_FSU   4 // self-unmapped
#define SAM_FMU   8 // mate-unmapped
#define SAM_FSR  16 // self on the reverse strand
#define SAM_FMR  32 // mate on the reverse strand
#define SAM_FR1  64 // this is read one
#define SAM_FR2 128 // this is read two
#define SAM_FSC 256 // secondary alignment


const char * read_pair[][5] = {
	{
		"read id",
		"TAACCCTAACCCTAACCCTAACCCTAACCCCAACCCCAACCCCAACCCCAACCCCAACCC",
		"BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB",
		"AAACCCTAACCCCTAACCCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAAC",
		"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"
	}
};

void check_error(int code, const char* error_msg)
{
	if (code)
	{
		fprintf(stderr, "%s\n", error_msg);
		fprintf(stderr, "\nError code: %d\n", code);
		abort();
	}
}

void print_alignments(const aln_batch* batch, FILE* out)
{
	for (int f = 0; f < batch->n_frags; ++f) {
		for (int i = 0; i < batch->n_reads_frag; ++i) {
			const aln_read*const read = aln_get_read(batch, f, i);

			fprintf(out, "%s", read->id);

			const aln_alignment*const aln = &read->alignment;
			int flag = 0;
			flag |= aln->paired         ? SAM_FPD : 0;
			flag |= aln->prop_paired    ? SAM_FPP : 0;
			flag |= aln->unmapped       ? SAM_FSU : 0;
			flag |= aln->reverse_strand ? SAM_FSR : 0;
			flag |= aln->secondary_aln  ? SAM_FSC : 0;

			fprintf(out, "\t%d", flag);
			if (aln->unmapped)
				fprintf(out, "\t*\t*");
			else
				fprintf(out, "\t%s\t%ld", aln->contig, aln->pos);

			fprintf(out, "\t%.*s", read->length, read->seq);
			if (read->qual)
				fprintf(out, "\t%.*s", read->length, read->qual);
			else
				fprintf(out, "\t*");
			// more stuff... insert size and tags
		}
	}
}


int main(int argc, const char* argv[])
{
	int error = 0;

	aln_opts opts;
	error = aln_init_opts(&opts);
	check_error(error, "Failed to init opts");

	//const aln_kv* tail = opts.parameters;
	//opts.parameters = aln_kv_insert_long("-t", 2, tail);
	//check_error(opts.parameters != NULL, "Failed to insert parameters");

	error = aln_init(&opts);
	check_error(error, "Failed to initialize");

	fprintf(stderr, "Initialized library version '%s'\n", aln_version());

	aln_ref ref;
	error = aln_load_ref("/u/pireddu/Projects/bwa_wrapper/mini_ref/mini_ref.fasta", &ref);
	check_error(error, "Failed to load reference");

	fprintf(stderr, "This reference has %d contigs\n", ref.n_contigs);

	aln_batch reads;
	error = aln_alloc_reads(&reads, 2, 1);//sizeof(read_pair) / (3 * sizeof(char*)));
	check_error(error, "Failed to allocate read batch");

	for (int f = 0; f < reads.n_frags; ++f) {
		const char* read_id = read_pair[f][0];
		for (int r = 0; r < reads.n_reads_frag; ++r) {
			fprintf(stderr, "setting read (%d, %d, %s, %s, %s\n", f, r, read_id, read_pair[f][1 + 2*r], read_pair[f][2 + 2*r]);
			error = aln_set_read(&reads, f, r, read_id, read_pair[f][1 + 2*r], read_pair[f][2 + 2*r], ALN_QUALITY_ENCODING_SANGER);
			check_error(error, "Failed to set read");
		}
	}

	fprintf(stderr, "Allocated read batch\n");
	fprintf(stderr, "\tn_frags: %d; n_reads_frag: %d\n",
			reads.n_frags, reads.n_reads_frag);

	//error = aln_align_reads(&reads, &opts);
	//check_error(error, "Failed to align reads!");

	//fprintf(stderr, "Reads aligned\n");

	aln_free_ref(&ref);
	aln_free_reads(&reads);
	aln_free_opts(&opts);

	return 0;
}
