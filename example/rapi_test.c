

#include <rapi.h>

#include <kstring.h>
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
		"GTCTCCCCCCAGGTGTGTGGTGATGCCAGGCATGCCCTTCCCCAGCATCAGGTCTCCAGAGCTGCAGAAGACGACGGCCGACTTGGATCACACTCTTGTG",
		"BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB",
		"GATGTCATCTGGAGCCCTGCTGCTTGCGGTGGCCTATAAAGCCTCCTAGTCTGGCTCCAAGGCCTGGCAGAGTCTTTCCCAGGGAAAGCTACAAGCAGCA",
		"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"
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

void print_alignments(const rapi_batch* batch, FILE* out)
{
	for (int f = 0; f < batch->n_frags; ++f) {
		for (int i = 0; i < batch->n_reads_frag; ++i) {
			const rapi_read*const read = rapi_get_read(batch, f, i);

			fprintf(out, "%s", read->id);

			for (int a = 0; a < read->n_alignments; ++a) {
				fprintf(out, "\tAlignment %d\t", a);
				const rapi_alignment*const aln = &read->alignments[a];
				int flag = 0;
				flag |= aln->paired         ? SAM_FPD : 0;
				flag |= aln->prop_paired    ? SAM_FPP : 0;
				flag |= (!aln->mapped)      ? SAM_FSU : 0;
				flag |= aln->reverse_strand ? SAM_FSR : 0;
				flag |= aln->secondary_aln  ? SAM_FSC : 0;

				fprintf(out, "\t%d", flag);
				if (aln->mapped)
					fprintf(out, "\t%s\t%ld", aln->contig->name, aln->pos);
				else
					fprintf(out, "\t*\t*");

				fprintf(out, "\t%.*s", read->length, read->seq);
				if (read->qual)
					fprintf(out, "\t%.*s", read->length, read->qual);
				else
					fprintf(out, "\t*");
				// more stuff... insert size and tags
			}
		}
	}
}

int main(int argc, const char* argv[])
{
	rapi_error_t error = 0;

	rapi_opts opts;
	error = rapi_opts_init(&opts);
	check_error(error, "Failed to init opts");

	//const rapi_kv* tail = opts.parameters;
	//opts.parameters = rapi_kv_insert_long("-t", 2, tail);
	//check_error(opts.parameters != NULL, "Failed to insert parameters");

	error = rapi_init(&opts);
	check_error(error, "Failed to initialize");

	fprintf(stderr, "Initialized API version '%s'\n", RAPI_API_VERSION);
	fprintf(stderr, "Using aligner '%s', version '%s'\n", rapi_aligner_name(), rapi_aligner_version());

	rapi_ref ref;
	error = rapi_ref_load("/u/pireddu/Projects/bwa_wrapper/mini_ref/mini_ref.fasta", &ref);
	check_error(error, "Failed to load reference");

	fprintf(stderr, "This reference has %d contigs\n", ref.n_contigs);

	rapi_batch reads;
	error = rapi_reads_alloc(&reads, 2, 1);//sizeof(read_pair) / (3 * sizeof(char*)));
	check_error(error, "Failed to allocate read batch");

	for (int f = 0; f < reads.n_frags; ++f) {
		const char* read_id = read_pair[f][0];
		for (int r = 0; r < reads.n_reads_frag; ++r) {
			fprintf(stderr, "setting read (%d, %d, %s, %s, %s)\n", f, r, read_id, read_pair[f][1 + 2*r], read_pair[f][2 + 2*r]);
			error = rapi_set_read(&reads, f, r, read_id, read_pair[f][1 + 2*r], read_pair[f][2 + 2*r], RAPI_QUALITY_ENCODING_SANGER);
			check_error(error, "Failed to set read");
		}
	}

	fprintf(stderr, "Allocated read batch\n");
	fprintf(stderr, "\tn_frags: %d; n_reads_frag: %d\n",
			reads.n_frags, reads.n_reads_frag);

	rapi_aligner_state* state;
	error = rapi_aligner_state_init(&opts, &state);
	check_error(error, "Failed to initialize aligner state");
	fprintf(stderr, "initialized aligner state\n");

	error = rapi_align_reads(&ref, &reads, 0, reads.n_frags, state);
	check_error(error, "Failed to align reads!");

	fprintf(stderr, "Reads aligned.  Now printing\n");

	kstring_t sam_buffer = { 0, 0, NULL };

	for (int frag = 0; frag < reads.n_frags; ++frag) {
		sam_buffer.l = 0;
		rapi_format_sam_b(&reads, frag, &sam_buffer);
		printf("%s\n", sam_buffer.s);
	}

	free(sam_buffer.s);
	rapi_aligner_state_free(state);
	rapi_ref_free(&ref);
	rapi_reads_free(&reads);
	rapi_opts_free(&opts);

	return 0;
}
