#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ntrx_expression_tools.h"
#include "ntrx_file_utils.h"
#include "ntrx_genomics_coords.h"
#include "ntrx_annot_region_tools.h"
#include "ntrx_human_data_processing.h"
#include "ntrx_gff_utils.h"
#include "ntrx_ansi_string.h"
#include <algorithm>
#include "ntrx_nomenclature.h"

int main(int argc, char* argv[])
{
	if(argc < 3)
	{
		fprintf(stderr, "USAGE: %s [options] [arguments]\n\
	-GENCODE_GTF_2_Interval_per_feature [GENCODE exons GTF file path] [Sub-element feature to merge (e.g., \"exon\")] [Super-element group string id to merge wrt (e.g., \"gene_id\")]\n\
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
	else if (strcmp(argv[1], "-GENCODE_GTF_2_Interval_per_feature") == 0)
	{
		//-GENCODE_GTF_2_Interval [GENCODE GTF file path]\n
		if (argc != 6)
		{
			fprintf(stderr, "USAGE: %s -GENCODE_GTF_2_Interval_per_feature [GENCODE exons GTF file path] [Sub-element feature to merge (e.g., \"exon\")] [Super-element group string id to merge wrt (e.g., \"gene_id\")] [Output file path]]\n", argv[0]);
			exit(0);
		}

		char* gtf_fp = argv[2];
		char* sub_element_feature = argv[3];
		char* super_element_grp_str_id = argv[4];
		char* op_fp = argv[5];

		// Load the gtf file.
		vector<t_annot_region*>* regions = load_GFF(gtf_fp);

		// Parse the group strings.
		parse_GENCODE_gff_grp_strs(regions);

		// Filter the regions per sub element feature type.
		vector<t_annot_region*>* sub_elements = get_gff_regions_per_feature_type(regions, sub_element_feature);
		if (sub_elements->size() == 0)
		{
			fprintf(stderr, "Could not find features of type %s, are you sure about the feature type?\n", sub_element_feature);
			exit(0);
		}

		vector<t_annot_region*>* merged_selected_regions = new vector<t_annot_region*>();
		vector<char*>* sel_region_props = new vector<char*>();

		// Form the largest regions.
		form_largest_regions_per_property_type(sub_elements,
			merged_selected_regions,
			sel_region_props,
			super_element_grp_str_id);

		fprintf(stderr, "Merged into %d super-regions.\n", (int)merged_selected_regions->size());

		vector<t_annot_region*>* merged_exonic_regions = new vector<t_annot_region*>();
		for (int i_merged_reg = 0; i_merged_reg < (int)merged_selected_regions->size(); i_merged_reg++)
		{
			// Chromosome check.
			vector<char*>* interval_chr_ids = get_chr_ids(merged_selected_regions->at(i_merged_reg)->intervals);
			if (interval_chr_ids->size() != 1)
			{
				fprintf(stderr, "** %s is located on multiple chromosomes. **\n", merged_selected_regions->at(i_merged_reg)->name);
			}

			// Merge all the CDSs for the current gene.
			vector<t_annot_region*>* merged_exons = merge_annot_regions(merged_selected_regions->at(i_merged_reg)->intervals, 0);

			// Get the introns for the current merged cds's.
			sort(merged_exons->begin(), merged_exons->end(), sort_regions);

			// Set the region.
			t_annot_region* cur_multi_exonic_region = get_empty_region();
			cur_multi_exonic_region->name = t_string::copy_me_str(sel_region_props->at(i_merged_reg));
			cur_multi_exonic_region->intervals = merged_exons;
			cur_multi_exonic_region->chrom = t_string::copy_me_str(merged_selected_regions->at(i_merged_reg)->chrom);
			cur_multi_exonic_region->start = merged_exons->at(0)->start;
			cur_multi_exonic_region->end = merged_exons->back()->end;
			cur_multi_exonic_region->strand = merged_selected_regions->at(i_merged_reg)->strand;
			cur_multi_exonic_region->intervals = merged_exons;

			// Add the current merged multiexonic region to the list of multiexonic regions.
			merged_exonic_regions->push_back(cur_multi_exonic_region);
			//delete_annot_regions(merged_exons);
		} // i_merged_reg loop.

		dump_Interval(op_fp, merged_exonic_regions);
	} // -GENCODE_GTF_2_Interval_per_feature option.
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

