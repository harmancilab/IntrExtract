#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ntrx_expression_tools.h"
#include "ntrx_file_utils.h"
#include "ntrx_genomics_coords.h"
#include "ntrx_annot_region_tools.h"
#include "ntrx_gff_utils.h"
#include "ntrx_ansi_string.h"
#include <algorithm>
#include "ntrx_nomenclature.h"

int main(int argc, char* argv[])
{
	if(argc < 3)
	{
		fprintf(stderr, "USAGE: %s [options] [arguments]\n\
	-compute_single_cell_expression_stats_per_10X_SAM\n\
	-get_sequence_mappability_stats_per_annotations\n", argv[0]);
		exit(0);
	}

	if (t_string::compare_strings(argv[1], "-get_sequence_mappability_stats_per_annotations"))
	{
		if (argc != 6)
		{
			fprintf(stderr, "USAGE: %s -get_sequence_mappability_stats_per_annotations [Annotation regions interval file path] [Multimapp dir.] [Genome sequence dir.] [Output file path]\n", argv[0]);
			exit(0);
		}

		char* annotation_regions_interval_fp = argv[2];
		char* multimapp_dir = argv[3];
		char* genome_seq_dir = argv[4];
		char* op_fp = argv[5];

		get_mappability_sequence_variates_per_intrexes(annotation_regions_interval_fp, multimapp_dir, genome_seq_dir, op_fp);
	} // -get_sequence_mappability_stats_per_annotations option.
	else if (t_string::compare_strings(argv[1], "-compute_single_cell_expression_stats_per_10X_SAM"))
	{
		if (argc != 6)
		{
			fprintf(stderr, "USAGE: %s -compute_single_cell_allelic_stats_per_10X_SAM [Cell barcodes list file path] [SAM file path (stdin ok)] [Annotation intervals file path] [Output file path]\n", argv[0]);
			exit(0);
		}

		char* per_cell_barcodes_fp = argv[2];
		char* SAM_fp = argv[3];
		char* annotation_regions_interval_fp = argv[4];
		char* op_fp = argv[5];

		compute_single_cell_expression_stats_per_10X_SAM(per_cell_barcodes_fp, SAM_fp, annotation_regions_interval_fp, op_fp);
	} // -compute_single_cell_expression_stats_per_10X_SAM option.
}

