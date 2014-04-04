

#include <aligner.h>

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


const char * read_pair[] = {
	"read id",
	"TAACCCTAACCCTAACCCTAACCCTAACCCCAACCCCAACCCCAACCCCAACCCCAACCC",
	"BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB",
	"TAACCCTAACCCCTAACCCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAAC",
	"BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB"
};

void check_error(int code, const char* error_msg) {
	if (code) {
		fprintf(stderr, error_msg);
		fprintf(stderr, "\nError code: %d\n", code);
		abort();
	}
}

void print_alignments(const aln_batch* batch, FILE* out) {
	for (int f = 0; f < batch->n_frags; ++f) {
		for (int i = 0; i < batch->n_reads_frags; ++i) {
			const aln_read* read = aln_get_read(batch, f, i);

			fprintf(out, read->id);

			const aln_alignment* aln = &read->alignment;
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

		}
	}
}


int main(int arc, const char* argv[]) {
	int error = 0;

	error = aln_init();
	check_error(error, "Failed to initialize");

	fprintf(stderr, "Initialized library version %s", aln_version());

	aln_opts opts;
	error = aln_init_opts(&opts);
	check_error(error, "Failed to init opts");

	const aln_kv* tail = opts.parameters;
	opts.parameters = aln_kv_insert_long("-t", 2, tail);

	check_error(opts.parameters != NULL, "Failed to insert parameters");

	aln_ref ref;
	error = aln_load_ref("/home/pireddu/Projects/bwa_wrapper/mini_ref/mini_ref", &ref);
	check_error(error, "Failed to load reference");


	aln_batch reads;
	error = aln_alloc_reads(&reads, 2, 1);//sizeof(read_pair) / (3 * sizeof(char*)));
	check_error(error, "Failed to allocate read batch");

	for (int f = 0; f < reads->n_frags; ++f) {
		const char* read_id = read_pair[f][0];
		for (int r = 1; r <= reads->n_reads_frags; ++r) {
			error = aln_set_read(f, r, read_id, read_pair[f][r], read_pair[f][r+1]);
			check_error(error, "Failed to set read");
		}
	}

	fprintf(stderr, "Allocated read batch\n");
	fprintf(stderr, "\tn_allocated_frags: %d; n_frags: %d; n_reads_frags: %d\n",
			reads.n_allocated_frags, reads.n_frags, reads.n_reads_frags);

	error = aln_align_reads(&reads, &opts);
	check_error(error, "Failed to align reads!");

	fprintf(stderr, "Reads aligned\n");

	check_error(error, "Failed to align reads!");
	aln_free_ref(&ref);
	aln_free_reads(&reads);

	return 0;
}
