#include <stdio.h>
#include <stdlib.h>
#include "ntrx_mapped_read_tools.h"
#include <vector>
#include <ctype.h>
#include <math.h>
#include "ntrx_signal_track_tools.h"
#include "ntrx_annot_region_tools.h"
#include "ntrx_genome_sequence_tools.h"
//#include "ntrx_variation_tools.h"
#include "ntrx_genomics_coords.h"
#include "ntrx_file_utils.h"
#include "ntrx_nomenclature.h"
#include "ntrx_nucleotide.h"
#include <string.h>
#include <algorithm>
#include "ntrx_ansi_string.h"

using namespace std; 

bool __DUMP_MAPPED_READ_TOOLS_MSGS__ = false;

FILE* get_processed_reads_ptr_wrapper(char* cur_chr_reads_fp)
{
	FILE* f_cur_chr_reads = open_f(cur_chr_reads_fp, "r");
//	if (t_string::compare_strings(cur_chr_reads_fp, "stdin"))
//	{
//		f_cur_chr_reads = stdin;
//	}
//	else if (t_string::ends_with(cur_chr_reads_fp, ".gz"))
//	{
//		char ungzip_cmd[1000];
//		sprintf(ungzip_cmd, "gzip -cd %s", cur_chr_reads_fp);
//#ifdef _WIN32
//		f_cur_chr_reads = _popen(ungzip_cmd, "r");
//#else 
//		f_cur_chr_reads = popen(ungzip_cmd, "r");
//#endif
//	}
//	else
//	{
//		f_cur_chr_reads = open_f(cur_chr_reads_fp, "r");
//	}

	return f_cur_chr_reads;
}


void count_10X_reads_per_cell_per_barcodes_list(char* TenX_SAM_file, char* per_cell_barcode_list_fp, char* op_fp)
{
	fprintf(stderr, "Loading and sorting barcodes from %s.\n", per_cell_barcode_list_fp);
	vector<char*>* per_cell_barcodes = buffer_file(per_cell_barcode_list_fp);
	sort(per_cell_barcodes->begin(), per_cell_barcodes->end(), t_string::sort_strings_per_prefix);
	int* n_reads_per_barcode = new int[per_cell_barcodes->size() + 2];
	memset(n_reads_per_barcode, 0, sizeof(int) * per_cell_barcodes->size());

	fprintf(stderr, "Generating the read count per cell for %d cell barcodes.\n", per_cell_barcodes->size());

	unsigned int n_reads_processed = 0;
	unsigned int n_unmatched_BC_reads = 0;
	FILE* f_sam = open_f(TenX_SAM_file, "r");
	while (1)
	{
		char* SAM_line = getline(f_sam);
		if (SAM_line == NULL)
		{
			break;
		}

		n_reads_processed++;
		if (n_reads_processed % (1000 * 1000) == 0)
		{
			fprintf(stderr, "Processing %d. read; %d unmatched BC reads.                       \r", n_reads_processed, n_unmatched_BC_reads);
		}

		// Parse the barcode and update counts.
		char barcode_buffer[100];
		if (get_10X_cell_barcode_per_SAM_read(SAM_line, barcode_buffer))
		{
			int cell_i = -1;
			int search_cell_i = t_string::fast_search_string_per_prefix(barcode_buffer, per_cell_barcodes, 0, per_cell_barcodes->size() - 1);
			while (search_cell_i > 0 &&
				(t_string::sort_strings_per_prefix(barcode_buffer, per_cell_barcodes->at(search_cell_i)) ||
					t_string::compare_strings(per_cell_barcodes->at(search_cell_i), barcode_buffer)))
			{
				search_cell_i--;
			} // search_cell_i loop.

			while (search_cell_i < per_cell_barcodes->size() &&
				(t_string::sort_strings_per_prefix(per_cell_barcodes->at(search_cell_i), barcode_buffer) ||
					t_string::compare_strings(per_cell_barcodes->at(search_cell_i), barcode_buffer)))
			{
				if (t_string::compare_strings(per_cell_barcodes->at(search_cell_i), barcode_buffer))
				{
					cell_i = search_cell_i;
					break;
				}
				else
				{
					search_cell_i++;
				}
			} // search_cell_i loop.

			  // If we found a cell, update the count of it.
			if (cell_i != -1)
			{
				n_reads_per_barcode[cell_i]++;
			}
			else
			{
				n_unmatched_BC_reads++;
			}
		} // tag parsing check.

		delete[] SAM_line;
	} // file reading loop.
	close_f(f_sam, TenX_SAM_file);

	// Write output.
	fprintf(stderr, "Saving the read counts per barcode to %s\n", op_fp);
	FILE* f_per_BC_read_counts = open_f(op_fp, "w");
	for (int bc_i = 0; bc_i < per_cell_barcodes->size(); bc_i++)
	{
		fprintf(f_per_BC_read_counts, "%s\t%d\n", per_cell_barcodes->at(bc_i), n_reads_per_barcode[bc_i]);
	} // bc_i loop.
	fclose(f_per_BC_read_counts);
}

bool get_10X_cell_barcode_per_SAM_read(char* TenX_SAM_read_line, char* barcode_buffer)
{
	char temp_tag_entry_buffer[1000];
	bool parsed_BC_tag = get_SAM_read_tag_entry(TenX_SAM_read_line, temp_tag_entry_buffer, "CB:Z:");

	// Copy the barcode.
	if (parsed_BC_tag)
	{
		strcpy(barcode_buffer, &(temp_tag_entry_buffer[5]));
	}

	return(parsed_BC_tag);
}


bool get_SAM_read_tag_entry(char* sam_read_line, char* tag_entry_buffer, const char* tag_prefix)
{
	int i_cur_char = 0;
	char cur_SAM_entry[1010];
	int l_entry_buff = 1000;

	// Find 11th token.
	int i_tok = 0;
	while (i_tok < 11)
	{
		t_string::get_next_token(sam_read_line, cur_SAM_entry, l_entry_buff, "\t", i_cur_char);
		i_tok++;
	} // i_tok loop.

	if (__DUMP_MAPPED_READ_TOOLS_MSGS__)
	{
		fprintf(stderr, "Located 11. token @ %d\n", i_cur_char);
	}

	// Loop through next set of tokens.
	bool found_tag = false;
	while (t_string::get_next_token(sam_read_line, cur_SAM_entry, l_entry_buff, "\t", i_cur_char))
	{
		// Copy the barcode sequence.
		if (t_string::starts_with(cur_SAM_entry, tag_prefix))
		{
			if (__DUMP_MAPPED_READ_TOOLS_MSGS__)
			{
				fprintf(stderr, "%s tag: %s\n", tag_prefix, cur_SAM_entry);
			}

			strcpy(tag_entry_buffer, cur_SAM_entry);
			found_tag = true;
			break;
		}
	} // i_tok loop. 

	return(found_tag);
}

void close_processed_reads_ptr_wrapper(FILE* f_cur_chr_reads, char* cur_chr_reads_fp)
{
	if (t_string::compare_strings(cur_chr_reads_fp, "stdin"))
	{
	}
	else if (t_string::ends_with(cur_chr_reads_fp, ".gz"))
	{
#ifdef _WIN32
		_pclose(f_cur_chr_reads);
#else 
		pclose(f_cur_chr_reads);
#endif
	}
	else
	{
		fclose(f_cur_chr_reads);
	}
}

bool set_preprocessed_read_file_path_per_dir_chr(char* preprocessed_reads_dir, char* chrom, char* preprocessed_reads_fp)
{
	sprintf(preprocessed_reads_fp, "%s/%s_mapped_reads.txt", preprocessed_reads_dir, chrom);

	if (check_file(preprocessed_reads_fp))
	{
		return true;
	}
	
	sprintf(preprocessed_reads_fp, "%s/%s_mapped_reads.txt.gz", preprocessed_reads_dir, chrom);

	if (check_file(preprocessed_reads_fp))
	{
		return true;
	}

	return false;
}

vector<int>* count_preprocessed_reads(char* preprocessed_reads_dir, vector<char*>* _chr_ids)
{
	// Load all the reads per chromosome: Load all the chromosomes.
	char chr_ids_list_fp[1000];
	sprintf(chr_ids_list_fp, "%s/chr_ids.txt", preprocessed_reads_dir);

	// Load the chromosome id's.
	vector<char*>* chr_ids = buffer_file(chr_ids_list_fp);
	if (chr_ids == NULL)
	{
		fprintf(stderr, "Could not load chromosome id's from %s\n", chr_ids_list_fp);
		return(NULL);
	}

	vector<int>* n_reads_per_chromosome = new vector<int>();

	int n_total_reads = 0;
	for (int i_chr = 0; i_chr < chr_ids->size(); i_chr++)
	{
		_chr_ids->push_back(t_string::copy_me_str(chr_ids->at(i_chr)));

		char cur_chr_reads_fp[1000];
		//sprintf(cur_chr_reads_fp, "%s/%s_mapped_reads.txt", preprocessed_reads_dir, chr_ids->at(i_chr));
		if (!set_preprocessed_read_file_path_per_dir_chr(preprocessed_reads_dir, chr_ids->at(i_chr), cur_chr_reads_fp))
		{
			fprintf(stderr, "Could not find preprocessed reads file for %s in %s.\n", preprocessed_reads_dir, chr_ids->at(i_chr));
			exit(0);
		}

		fprintf(stderr, "Counting the reads in %s (%s)\n", chr_ids->at(i_chr), cur_chr_reads_fp);

		//int n_cur_chr_reads = 0;
		FILE* f_cur_chr_reads = get_processed_reads_ptr_wrapper(cur_chr_reads_fp);

		int n_cur_chr_reads = get_n_lines(f_cur_chr_reads);

		close_processed_reads_ptr_wrapper(f_cur_chr_reads, cur_chr_reads_fp);

		n_reads_per_chromosome->push_back(n_cur_chr_reads);

		// Update total # of reads.
		n_total_reads += n_cur_chr_reads;

		//fprintf(stderr, "%s: %d reads. (%d reads in total)\n", chr_ids->at(i_chr), n_cur_chr_reads, n_total_reads);
	} // i_chr loop.

	return(n_reads_per_chromosome);
}

bool sort_read_line_entries_per_id(t_read_line_w_id* read1, t_read_line_w_id* read2)
{
	return(t_string::sort_strings(read1->id, read2->id));
}

void preprocess_SAM_read_line(char* cur_line,
	char* read_id,
	char* chrom,
	int& chr_index, int& sequenced_length,
	char& strand_char,
	char* cigar_str)
{
	// Skip the comment and headers.
	if (cur_line[0] == '@')
	{
		chrom[0] = 0;
		chr_index = 0;
		return;
	}

	int flag;
	char flag_str[100];
	int _chr_index;
	char _chr_index_str[100];
	char fragment[100000];
	char phred_quality_str[100000];

	//t_string_tokens* cur_tokens = t_string::tokenize_by_chars(cur_line, "\t");
	//if(sscanf(cur_line, "%s %d %s %d %*s %s %*s %*s %*s %s %s", read_id, &flag, chrom, &_chr_index, cigar_str, fragment, phred_quality_str) == 7)
	//if(cur_tokens->size() >= 11)
	if (sscanf(cur_line, "%[^\t] %[^\t] %[^\t] %[^\t] %*[^\t] %[^\t] %*[^\t] %*[^\t] %*[^\t] %[^\t] %[^\t]", read_id, flag_str, chrom, _chr_index_str, cigar_str, fragment, phred_quality_str) == 7)
	{
		//t_string::copy(read_id, cur_tokens->at(0)->str());
		//flag = atoi(cur_tokens->at(1)->str());
		//t_string::copy(chrom, cur_tokens->at(2)->str());
		//_chr_index = atoi(cur_tokens->at(3)->str());
		//t_string::copy(cigar_str, cur_tokens->at(5)->str());
		//t_string::copy(fragment, cur_tokens->at(9)->str());
		//t_string::copy(phred_quality_str, cur_tokens->at(10)->str());

		_chr_index = atoi(_chr_index_str);
		flag = atoi(flag_str);

		// Translate the 0 based index in SAM file to ELAND's 1 based indexing.
		//chr_index += (CODEBASE_START_BASE - SAM_START_BASE);
		_chr_index = translate_coord(_chr_index, SAM_COORDS::start_base, CODEBASE_COORDS::start_base);

		// Check the flag and determine the strand.
		strand_char = 'F';
		if (flag & 0x10)
		{
			strand_char = 'R';
		}

		// Sanity check. Is this fragment mapped?
		if (flag & 0x04)
		{
			// The read is not mapping.
			chrom[0] = 0;
		}
		else
		{
			chr_index = _chr_index;
			sequenced_length = strlen(fragment);
		}
	}
	else
	{
		chrom[0] = 0;
	}

	//t_string::clean_tokens(cur_tokens);
}


#define __UCHAR_MAPPABILITY__

double* load_normalized_multimappability_profile(char* mapability_signal_profile_fp, int& l_mapability_profile)
{
	double* mapability_signal_profile = NULL;

	if(!check_file(mapability_signal_profile_fp))
	{
		l_mapability_profile = 0;
		return(NULL);
	}

#ifdef __DOUBLE_MAPPABILITY__
	// Load the mapability map signal profile, do filtering.
	mapability_signal_profile = load_per_nucleotide_binary_profile(mapability_signal_profile_fp, l_mapability_profile);

	// Do mapability aware median filtering on the current signal profile.
if(__DUMP_PEAK_CALLING_UTILS_MSGS__)
	fprintf(stderr, "Scaling the mapability map with %d.\n", l_read_mapability_signal * 2);

	int mapability_scaling = l_read * 2;
	for(int i = 1; i <= l_mapability_profile; i++)
	{
		mapability_signal_profile[i] /= mapability_scaling;
	} // i loop.
#elif defined(__UCHAR_MAPPABILITY__)
	// Following loads the mappability signal profile from the char version of the multi-mappability profile.
	// Load the mapability map signal profile, do filtering.
	unsigned char* mapability_signal_char_profile = load_per_nucleotide_binary_uchar_profile(mapability_signal_profile_fp, l_mapability_profile);
	mapability_signal_profile = new double[l_mapability_profile + 2];
	for(int i = 1; i <= l_mapability_profile; i++)
	{
		unsigned char unsigned_char_val = (unsigned char)(mapability_signal_char_profile[i]);
		mapability_signal_profile[i] = (double)(unsigned_char_val);
		mapability_signal_profile[i] /= 100;

		if(mapability_signal_profile[i] < 0)
		{
			fprintf(stderr, "Sanity check failed.\n");
			exit(0);
		}
	} // i loop.
	delete [] mapability_signal_char_profile;
#else
	#error "Must define the type of mappability."
#endif

	return(mapability_signal_profile);
}

bool sort_read_lines(char* read1, char* read2)
{
	return(t_string::sort_strings(read1, read2));
}

#define MAX(x, y) (((x)>(y))?(x):(y))
#define MIN(x, y) (((x)<(y))?(x):(y))

int get_chr_i_per_group_entry(int self_flag_w_chr_w_mate_chr_val, bool self)
{
	if(self)
	{
		int flagged_val = (self_flag_w_chr_w_mate_chr_val & 0xFF00);
		int self_chr_i = flagged_val >> 8;
		return(self_chr_i);
	}
	else
	{
		int flagged_val = (self_flag_w_chr_w_mate_chr_val & 0xFF);
		int mate_chr_i = flagged_val;
		return(mate_chr_i);
	}
}

unsigned int get_flag_per_group_entry(int self_flag_w_chr_w_mate_chr_val)
{
	int flagged_val = (self_flag_w_chr_w_mate_chr_val & 0xFFFF0000);
	unsigned int flag_val = flagged_val >> 16;
	return(flag_val);
}

bool sort_barcode_entries_per_barcode(t_10X_BX_barcode_group_entry* ent1, t_10X_BX_barcode_group_entry* ent2)
{
	return(ent1->BX_barcode < ent2->BX_barcode);
}

int recursive_locate_barcode_entry_index(double cur_barcode, 
										vector<t_10X_BX_barcode_group_entry*>* barcode_entries, 
										int start_i, int end_i)
{
	int mid_i = (start_i + end_i) / 2;

	if(mid_i == start_i || mid_i == end_i)
	{
		return(mid_i);
	}

	if(cur_barcode > barcode_entries->at(mid_i)->BX_barcode)
	{
		return recursive_locate_barcode_entry_index(cur_barcode, 
													barcode_entries, 
													mid_i, end_i);
	}
	else if(cur_barcode < barcode_entries->at(mid_i)->BX_barcode)
	{
		return recursive_locate_barcode_entry_index(cur_barcode, 
													barcode_entries, 
													start_i, mid_i);
	}
	else
	{
		return mid_i;
	}

}

double get_BX_barcode_val_per_barcode_str(char* barcode_str, int* per_char_val, double base_val, int& n_processed_chars)
{
	double cur_barcode = 0;
	double cur_base = 1;
	int i = 0;
	while(barcode_str[i] &&
		per_char_val[(int)(barcode_str[i])] < 4)
	{
		int num_val = per_char_val[(int)(barcode_str[i])];
		cur_barcode += (cur_base * num_val);

		cur_base *= base_val;

		i++;
	} // i loop.

	n_processed_chars = i;

	return(cur_barcode);
}


void get_barcode_string_per_value(double barcode_val, double base_val, char* barcode_str)
{
	char nuc_per_val[] = "ACGTN";

	double cur_val = barcode_val;
	int nuc_i = 0;
	while(cur_val > 0)
	{
		int cur_nuc_val = cur_val - floor(cur_val / base_val) * base_val;
		barcode_str[nuc_i] = nuc_per_val[cur_nuc_val];
		cur_val = floor(cur_val / base_val);

		nuc_i++;
	} // cur_val loop.
	barcode_str[nuc_i] = 0;

	//fprintf(stderr, "Uninverted: %s\n", barcode_str);
	//t_string::revert(barcode_str);
}

void get_unique_10X_BX_tags(char* tenX_sam_fp, int min_mapQ, char* op_fp)
{
	fprintf(stderr, "Parsing BX tags from 10X reads in %s using minimum mapping quality of %d\n", tenX_sam_fp, min_mapQ);

	FILE* f_tenX_sam = NULL;

	if(t_string::compare_strings(tenX_sam_fp, "stdin"))
	{
		f_tenX_sam = stdin;
	}
	else
	{
		f_tenX_sam = open_f(tenX_sam_fp, "r");
	}

	int* per_char_val = get_per_char_as_nuc_2_num_coding_array();

	int n_reads_no_BX = 0;
	int n_processed_reads = 0;

	int l_barcode = 16;

	// This is the list of barcode groups.
	vector<double>* bx_tag_vals = new vector<double>();
	while(1)
	{
		char* cur_line = getline(f_tenX_sam);
		if(cur_line == NULL)
		{
			break;
		}

		n_processed_reads++;

		if(n_processed_reads % (1000*1000) == 0)
		{
			fprintf(stderr, "Processing %d. read. (%d w BX, %d w/o BX)                    \r", n_processed_reads, (int)bx_tag_vals->size(), n_reads_no_BX);
		}

		int flag = 0;
		char read_id[1000];
		char flag_str[100];
		int _chr_index = 0;
		char _chr_index_str[100];
		int _mate_chr_index = 0;
		char _mate_chr_index_str[100];
		char _mapQ_str[100];
		char fragment[100000];
		char phred_quality_str[100000];
		char cigar_str[1000];
		char chrom[100];
		char mate_chrom[100];

		// Parse the SAM line.
		if(sscanf(cur_line, "%[^\t] %[^\t] %[^\t] %[^\t] %[^\t] %[^\t] %[^\t] %[^\t] %*[^\t] %[^\t] %[^\t]", 
							read_id, flag_str, chrom, _chr_index_str, _mapQ_str, cigar_str, 
							mate_chrom, _mate_chr_index_str, 
							fragment, phred_quality_str) == 10)
		{
			int mapQ = atoi(_mapQ_str);

			if(mapQ < min_mapQ)
			{

			}
			else
			{
				// Parse this read.
				bool found_BX = false;

				char cur_entry[1000];
				int i_cur_char = 0;
				while(t_string::get_next_token(cur_line, cur_entry, 1000, "\t", i_cur_char))
				{
					if(cur_entry[0] == 'B' && t_string::starts_with(cur_entry, "BX:Z:"))
					{			
						// Create and search the current list of GEM groups.
						char* cur_full_barcode_str = cur_entry;
						int n_proc_chars = 0;
						double cur_barcode = get_BX_barcode_val_per_barcode_str(&(cur_full_barcode_str[5]), per_char_val, 4, n_proc_chars);
						if(n_proc_chars != l_barcode)
						{
							fprintf(stderr, "Could not read %d characters from barcode: %s\n", l_barcode, &(cur_full_barcode_str[5]));
						}

						//char decoded_str[1000];
						//get_barcode_string_per_value(cur_barcode, 4, decoded_str);

						//if(!t_string::compare_strings(decoded_str, &(cur_full_barcode_str[5])))
						//{
						//	fprintf(stderr, "Decoded does not match original:\n%s\n%s\n", 
						//			decoded_str, &(cur_full_barcode_str[5]));
						//	exit(0);
						//}

						// Find this barcode among existing barcodes.
						bx_tag_vals->push_back(cur_barcode);
						found_BX = true;
						break;
					}
				} // token loop.

				if(!found_BX)
				{
					//fprintf(stderr, "Could not parse GEM group entry: %s\n", cur_line);
					n_reads_no_BX++;
				}
			} // mapQ check.
		} // Read info parse check.

		// Free memory.
		delete [] cur_line;
	} // file reading loop.

	if(t_string::compare_strings(tenX_sam_fp, "stdin"))
	{
	}
	else
	{
		fclose(f_tenX_sam);
	}

	fprintf(stderr, "%d reads with BX tags (%d reads missing BX), saving to %s\n", (int)bx_tag_vals->size(), n_reads_no_BX, op_fp);

	// Dump the unique tag values.
	sort(bx_tag_vals->begin(), bx_tag_vals->end());

	FILE* f_op = open_f(op_fp, "w");
	double cur_tag_val = -1;
	for(int i_tag = 0; i_tag < (int)bx_tag_vals->size(); i_tag++)
	{
		if(cur_tag_val != bx_tag_vals->at(i_tag))
		{
			fprintf(f_op, "%.1f\n", bx_tag_vals->at(i_tag));			
		}

		cur_tag_val = bx_tag_vals->at(i_tag);
	} // i_tag loop.
	fclose(f_op);
}

// Parse the BX tag::https://support.10xgenomics.com/genome-exome/software/pipelines/latest/output/bam
void parse_10x_linked_reads_per_BX_tag(char* tenX_sam_fp, vector<char*>* chr_ids, char* GEM_barcodes_list_fp, int min_mapQ, char* op_fp)
{
	fprintf(stderr, "Parsing linked reads per BX tag from 10X sam file @ %s (min mapQ %d) using %d chromosomes with barcodes in file %s.\n", tenX_sam_fp, min_mapQ, (int)chr_ids->size(), GEM_barcodes_list_fp);

	FILE* f_tenX_sam = NULL;

	if(t_string::compare_strings(tenX_sam_fp, "stdin"))
	{
		f_tenX_sam = stdin;
	}
	else
	{
		f_tenX_sam = open_f(tenX_sam_fp, "r");
	}

	int* per_char_val = get_per_char_as_nuc_2_num_coding_array();

	// Load the list of barcode groups.
	vector<t_10X_BX_barcode_group_entry*>* barcode_entries = new vector<t_10X_BX_barcode_group_entry*>();
	fprintf(stderr, "Loading GEM/BX barcodes.\n");
	vector<char*>* barcodes = buffer_file(GEM_barcodes_list_fp);
	if(barcodes == NULL)
	{
		fprintf(stderr, "Could not load bardcodes from %s\n", GEM_barcodes_list_fp);
		exit(0);
	}

	for(int ibc = 0; ibc < (int)barcodes->size(); ibc++)
	{
		t_10X_BX_barcode_group_entry* new_entry = new t_10X_BX_barcode_group_entry();
		new_entry->BX_barcode = atof(barcodes->at(ibc));
		new_entry->self_chr_w_mate_chr_w_flag = new vector<int>();
		new_entry->posn = new vector<int>();
		new_entry->mate_posn = new vector<int>();

		barcode_entries->push_back(new_entry);
	} // ibc loop.
	fprintf(stderr, "Loaded %d entries.\n", (int)barcode_entries->size());

	// Sort the entries.
	sort(barcode_entries->begin(), barcode_entries->end(), sort_barcode_entries_per_barcode);

	int l_barcode = 16;

	int n_processed_reads = 0;
	fprintf(stderr, "Loading and parsing 10X SAM file.\n");
	while(1)
	{
		char* cur_line = getline(f_tenX_sam);
		if(cur_line == NULL)
		{
			break;
		}

		n_processed_reads++;

		if(n_processed_reads % 1000000 == 0)
		{
			fprintf(stderr, "Processing %d. read: %d barcode entries                  \r", n_processed_reads, (int)barcode_entries->size());
		}

		// Parse this read.
		int flag = 0;
		char read_id[1000];
		char flag_str[100];
		int _chr_index = 0;
		char _chr_index_str[100];
		int _mate_chr_index = 0;
		char _mate_chr_index_str[100];
		char _mapQ_str[100];
		char fragment[100000];
		char phred_quality_str[100000];
		char cigar_str[1000];
		char chrom[100];
		char mate_chrom[100];

		unsigned char cur_chrom_i = (int)chr_ids->size();
		unsigned char mate_chrom_i = (int)chr_ids->size();

		// Parse the SAM line.
		if(sscanf(cur_line, "%[^\t] %[^\t] %[^\t] %[^\t] %[^\t] %[^\t] %[^\t] %[^\t] %*[^\t] %[^\t] %[^\t]", 
							read_id, flag_str, chrom, _chr_index_str, _mapQ_str, cigar_str, 
							mate_chrom, _mate_chr_index_str, 
							fragment, phred_quality_str) == 10)
		{
			int mapQ = atoi(_mapQ_str);

			if(mapQ >= min_mapQ)
			{
				_chr_index = atoi(_chr_index_str);
				_mate_chr_index = atoi(_mate_chr_index_str);
				flag = atoi(flag_str);		

				// Translate the SAM indexing.
				_chr_index = translate_coord(_chr_index, SAM_COORDS::start_base, CODEBASE_COORDS::start_base);
				_mate_chr_index = translate_coord(_mate_chr_index, SAM_COORDS::start_base, CODEBASE_COORDS::start_base);

				// Check the flag and determine the strand.
				char strand_char = 'F';
				if(flag & 0x10)
				{
					strand_char = 'R';
				}

				normalize_chr_id(chrom);							

				if(mate_chrom[0] == '=')
				{
					strcpy(mate_chrom, chrom);
				}
				else
				{
					normalize_chr_id(mate_chrom);
				}

				unsigned char mate_chrom_i = (int)chr_ids->size();
				mate_chrom_i = t_string::get_i_str(chr_ids, mate_chrom);

				// Sanity check. Is this fragment mapped?
				unsigned char cur_chrom_i = t_string::get_i_str(chr_ids, chrom);

				// Make sure this read passes the quality checks.
				if((flag & 0x400) == 0 &&					// Optical OR PCR duplicate.
					(flag & 0x200) == 0 &&					// Vendor quality checks.
					(flag & 0x4) == 0 &&					// Read mapped
					cur_chrom_i != (int)chr_ids->size())
				{
					// Find the group that this read belongs to.
					bool found_BX = false;

					// Find the "BX:Z:" substring index.
					int char_i = 0;
					while(cur_line[char_i])
					{
						if(cur_line[char_i] == 'B' &&
							cur_line[char_i+1] == 'X' &&
							cur_line[char_i+2] == ':' &&
							cur_line[char_i+3] == 'Z' &&
							cur_line[char_i+4] == ':')
						{
							found_BX = true;

							//// Find the next tab and set it to 0.
							//for(int char_i2 = char_i+5; 
							//	cur_line[char_i2] != 0; 
							//	char_i2++)
							//{
							//	if(cur_line[char_i2] == '\t')
							//	{
							//		cur_line[char_i2] = 0;
							//	}
							//} // char_i2 loop.
							
							break;
						} // BX:Z: check.
						char_i++;
					} // char_i loop.

					if(found_BX)
					{
						// Create and search the current list of GEM groups.
						// Following fixes the naming for the barcode.
						char* cur_full_barcode_str = &(cur_line[char_i]);
						char* cur_barcode_str = &(cur_full_barcode_str[5]);

						int n_proc_chars = 0;
						double cur_barcode = get_BX_barcode_val_per_barcode_str(cur_barcode_str, per_char_val, 4, n_proc_chars);
						if(n_proc_chars != l_barcode)
						{
							fprintf(stderr, "Could not read %d (%d) characters from barcode: %s\n", l_barcode, n_proc_chars, cur_barcode_str);
							exit(0);
						}

						// Find this barcode among existing barcodes.
						int ent_i = recursive_locate_barcode_entry_index(cur_barcode, barcode_entries, 0, (int)barcode_entries->size() - 1);

						while(ent_i > 0 && 
							ent_i < (int)barcode_entries->size() &&
							barcode_entries->at(ent_i)->BX_barcode >= cur_barcode)
						{
							ent_i--;
						} // ent_i rewind loop.

						bool found_matching_barcode_entry = false;
						while(ent_i >= 0 && 
							ent_i < (int)barcode_entries->size() &&
							barcode_entries->at(ent_i)->BX_barcode <= cur_barcode)
						{
							if(barcode_entries->at(ent_i)->BX_barcode == cur_barcode)
							{
								found_matching_barcode_entry = true;

								// Pack some values.
								int self_flag_w_chr_w_mate_chr_val = 0;
								self_flag_w_chr_w_mate_chr_val = (flag << 16) | (cur_chrom_i << 8) | (mate_chrom_i);

								barcode_entries->at(ent_i)->self_chr_w_mate_chr_w_flag->push_back(self_flag_w_chr_w_mate_chr_val);
								barcode_entries->at(ent_i)->posn->push_back(_chr_index);
								barcode_entries->at(ent_i)->mate_posn->push_back(_mate_chr_index);

								// Test if we can successfully retrieve the values.
								int rec_flag = get_flag_per_group_entry(self_flag_w_chr_w_mate_chr_val);
								int rec_self_chr_i = get_chr_i_per_group_entry(self_flag_w_chr_w_mate_chr_val, true);
								int rec_mate_chr_i = get_chr_i_per_group_entry(self_flag_w_chr_w_mate_chr_val, false);
								if(rec_flag != flag)
								{
									fprintf(stderr, "Recovered flag does not match: %d, %d: %s\n", rec_flag, (int)flag, cur_line);
									exit(0);
								}

								if(rec_self_chr_i != cur_chrom_i ||
									rec_mate_chr_i != mate_chrom_i)
								{
									fprintf(stderr, "Recovered chromosome indices do not match: %d, %d; %d, %d\n%s\n", 
										(int)rec_self_chr_i, (int)cur_chrom_i, 
										(int)rec_mate_chr_i, (int)mate_chrom_i,
										cur_line);
									exit(0);
								}

								// Break the barcode search loop.
								break;
							} // barcode match check.

							ent_i++;
						} // ent_i rewind loop.

						if(found_matching_barcode_entry == false)
						{
							fprintf(stderr, "Could not find entry with barcode %lf\n", cur_barcode);
						}
					} // BX:Z check.

					if(!found_BX)
					{
						//fprintf(stderr, "Could not parse GEM group entry: %s\n", cur_line);
					}
				} // read unmapping check from flag.
			} // Minimum mapQ check.
		} // sam file parse check.
		else
		{
			fprintf(stderr, "Could not parse %s\n", cur_line);
			exit(0);
		}

		delete [] cur_line;
	} // file reading loop.

	// Save the linked read information.	

	// Dump the text output.
#define __DUMP_TEXT__
#undef __DUMP_BINARY__

	FILE* f_op = NULL;
#ifdef __DUMP_TEXT__
	fprintf(stderr, "Saving text linked read information to %s\n", op_fp);
	f_op = open_f(op_fp, "w");
#elif defined(__DUMP_BINARY__)
	fprintf(stderr, "Saving binary linked read information to %s\n", op_fp);
	f_op = open_f(op_fp, "wb");
#endif 

	for(int bc_i = 0; bc_i < (int)barcode_entries->size(); bc_i++)
	{
		if((int)barcode_entries->at(bc_i)->self_chr_w_mate_chr_w_flag->size() <= 1)
		{
			continue;
		}

#ifdef __DUMP_TEXT__
		fprintf(f_op, "%.1f\t%d", barcode_entries->at(bc_i)->BX_barcode, (int)barcode_entries->at(bc_i)->self_chr_w_mate_chr_w_flag->size());
#elif defined(__DUMP_BINARY__)
		// Write the barcode.
		fwrite(&(barcode_entries->at(bc_i)->BX_barcode), sizeof(double), 1, f_op);

		// Write the number of entries in the barcode group.
		int n_reads = (int)(barcode_entries->at(bc_i)->self_chr_w_mate_chr_w_flag->size());
		fwrite(&(n_reads), sizeof(int), 1, f_op);
#endif 
		for(int i_r = 0; 
			i_r < (int)barcode_entries->at(bc_i)->self_chr_w_mate_chr_w_flag->size();
			i_r++)
		{
#ifdef __DUMP_TEXT__
			fprintf(f_op, "\t%d\t%d\t%d", 
					barcode_entries->at(bc_i)->self_chr_w_mate_chr_w_flag->at(i_r), 
					barcode_entries->at(bc_i)->posn->at(i_r), 
					barcode_entries->at(bc_i)->mate_posn->at(i_r));
#elif defined(__DUMP_BINARY__)
			int cur_self_chr_w_mate_chr_w_flag = barcode_entries->at(bc_i)->self_chr_w_mate_chr_w_flag->at(i_r);
			int cur_posn = barcode_entries->at(bc_i)->posn->at(i_r);
			int cur_mate_posn = barcode_entries->at(bc_i)->mate_posn->at(i_r);

			// Write the information about this entry.
			fwrite(&(cur_self_chr_w_mate_chr_w_flag), sizeof(int), 1, f_op);
			fwrite(&(cur_posn), sizeof(int), 1, f_op);
			fwrite(&(cur_mate_posn), sizeof(int), 1, f_op);
#endif 
		} // i_r loop.

#ifdef __DUMP_TEXT__
		fprintf(f_op, "\n");
#endif
	} // bc_i loop.
	fclose(f_op);

	if(t_string::compare_strings(tenX_sam_fp, "stdin"))
	{
	}
	else
	{
		fclose(f_tenX_sam);
	}
} // parse_10x_linked_reads_per_BX_tag

vector<t_annot_region*>* load_per_chrom_per_bin_discordancy_stats(char* per_bin_discordancy_fp, 
																					vector<char*>* chr_ids,
																					int l_bin)
{
	vector<t_annot_region*>* per_bin_discordancy_stat_regs = new vector<t_annot_region*>();

	FILE* f_per_bin_discordancy = open_f(per_bin_discordancy_fp, "r");
	while(1)
	{
		char* cur_line = getline(f_per_bin_discordancy);
		if(cur_line == NULL)
		{
			break;
		}

		char cur_chr_id[1000];
		int cur_bin_start = 0;
		double avg_multimapp = 0;
		int n_cis_discordant = 0;
		int n_trans_discordant = 0;
		if(sscanf(cur_line, "%s %d %lf %d %d", cur_chr_id, &cur_bin_start, &avg_multimapp, &n_cis_discordant, &n_trans_discordant) != 5)
		{
			fprintf(stderr, "Could not parse.\n");
			exit(0);
		}

		t_per_bin_10x_discordancy_stat* discordant_stats = new t_per_bin_10x_discordancy_stat();
		discordant_stats->cur_bin_start = cur_bin_start;		
		discordant_stats->avg_multimapp = avg_multimapp;
		discordant_stats->n_cis_discordant_reads = n_cis_discordant;
		discordant_stats->n_trans_discordant_reads = n_trans_discordant;

		t_annot_region* cur_reg = get_empty_region();
		cur_reg->chrom = t_string::copy_me_str(cur_chr_id);

		// We put a padding to get around intersecting of neighboring regions.
		cur_reg->start = cur_bin_start + 10;
		cur_reg->end = cur_bin_start + l_bin - 1 - 10;
		cur_reg->strand = '+';
		cur_reg->data = discordant_stats;
		per_bin_discordancy_stat_regs->push_back(cur_reg);
	} // file reading loop.
	fclose(f_per_bin_discordancy);

	return(per_bin_discordancy_stat_regs);
}

vector<int>* get_chromosome_lengths_per_mapped_reads(char* mapped_reads_dir)
{
	vector<int>* chr_lengths = new vector<int>();

	char chr_ids_fp[1000];
	sprintf(chr_ids_fp, "%s/chr_ids.txt", mapped_reads_dir);
	vector<char*>* chr_ids = buffer_file(chr_ids_fp);
	for(int i_chr = 0; i_chr < (int)chr_ids->size(); i_chr++)
	{
		int l_cur_chr = 0;
		char cur_line[1000];
		char cur_mapped_reads_fp[1000];
		sprintf(cur_mapped_reads_fp, "%s/%s_mapped_reads.txt", mapped_reads_dir, chr_ids->at(i_chr));
		FILE* f_mapped_reads = get_processed_reads_ptr_wrapper(cur_mapped_reads_fp);

		while(1)
		{
			if(fgets(cur_line, 1000, f_mapped_reads) == NULL)
			{
				break;
			}

			int cur_pos = 0;
			sscanf(cur_line, "%*s %*s %d", &cur_pos);

			if(cur_pos > l_cur_chr)
			{
				l_cur_chr = cur_pos + 1000;
			}
		} // file reading loop.
		//fclose(f_mapped_reads);
		close_processed_reads_ptr_wrapper(f_mapped_reads, cur_mapped_reads_fp);

		chr_lengths->push_back(l_cur_chr);
	} // i_chr loop.

	return(chr_lengths);
}

// Following is for sorting the mapped reads offline.
bool sort_read_line_info(t_read_line_sorting_info* info1, t_read_line_sorting_info* info2)
{
	return(info1->start < info2->start);
}

bool sort_regions_coords_first_names_second(t_annot_region* reg1, t_annot_region* reg2)
{
	if (reg1->start != reg2->start)
	{
		return(reg1->start < reg2->start);
	}
	else if (reg1->end != reg2->end)
	{
		return(reg1->end < reg2->end);
	}
	else
	{
		// Both starts and ends are the same; check the alleles.
		return t_string::sort_strings(reg1->name, reg2->name);
	}
}

enum{VAL, TYPE};
bool validate_mapping_map_str(char* mapping_map_str, bool& is_read_spliced)
{
	int i = 0;

	is_read_spliced = false;

	int state = VAL;
	while(mapping_map_str[i] != 0)
	{		
		if(state == VAL)
		{
			//fprintf(stderr, "%c (%d)\n", quality_str[i], 0);
			// MIDNSHPX=
			if(mapping_map_str[i] == 'M' ||
				mapping_map_str[i] == 'I' ||
				mapping_map_str[i] == 'D' ||
				mapping_map_str[i] == 'N' ||
				mapping_map_str[i] == 'S' ||		
				mapping_map_str[i] == 'H' ||
				mapping_map_str[i] == 'P' ||
				mapping_map_str[i] == 'X' ||
				mapping_map_str[i] == '=')
			{
				state = TYPE;

				//if(mapping_map_str[i] != 'M')
				// If the state is N, we assume that this is a spliced read: Page 7, CIGAR definition @ https://samtools.github.io/hts-specs/SAMv1.pdf 
				if (mapping_map_str[i] == 'N')
				{
					is_read_spliced = true;
				}
			}
			else if(mapping_map_str[i] >= '0' && mapping_map_str[i] <= '9')
			{
				// State is still VAL.
			}
			else
			{
				return(false);
			}
		}
		else if(state == TYPE)
		{
			//fprintf(stderr, "%c (%d)\n", quality_str[i], 1);
			// A number is expected.
			if(mapping_map_str[i] >= '0' && mapping_map_str[i] <= '9')
			{
				state = VAL;
			}
			else
			{
				return(false);
			}
		}

		// Move to next character.
		i++;
	}

	return(true);
}

void get_next_entry_per_mapp_map_string(char* mapping_map_str,
										int& i_mapp_map, 
										bool& is_matching,
										//t_string* cur_entry_length_str,
										int& l_cur_entry,
										char& entry_type_char)
{	
	// Clean the length string.
	//cur_entry_length_str->empty();
	l_cur_entry = 0;

	// Get the next entry in the cigar string.
	while(mapping_map_str[i_mapp_map] != 0)
	{
		if(mapping_map_str[i_mapp_map] < '0' || mapping_map_str[i_mapp_map] > '9')
		{
			break;
		}
		//cur_entry_length_str->concat_char(mapping_map_str[i_mapp_map]);
		l_cur_entry = l_cur_entry*10 + (int)(mapping_map_str[i_mapp_map]-'0');
		i_mapp_map++;
	}

	is_matching = false;
	if(mapping_map_str[i_mapp_map] == 'M')
	{
		//fprintf(stderr, "Adding matching length of %d\n", l_cur_entry);
		is_matching = true;
	}
	else
	{
		//fprintf(stderr, "Adding some other length of %d\n", l_cur_entry);
	}	

	entry_type_char = mapping_map_str[i_mapp_map];

	// Move over the current entry identifier.
	i_mapp_map++;
}

// This is a generic iterator over the processed reads file, useful for doing online processing of the reads files, for example, when they are too large to load into memory.
void preprocessed_read_file_iterator(char* mapped_reads_fp,
	void (per_read_callback)(char*, char, int, void*), 
	void (per_fragment_callback)(char*, char, int, void*),
	void* per_read_callback_param,
	void* per_fragment_callback_param)
{
	printf("Loading mapped-reads from %s.\n", mapped_reads_fp);

	FILE* f_mrf = open_f(mapped_reads_fp, "r");

	//char cur_fragment[10000];
	char mapping_map_str[10000];
	char strand_char;
	int chr_index;

	// Read and validate the mapped reads in the file.
	//while(fscanf(f_mapped_reads, "%s %s %c %d", cur_fragment, quality_str, &strand_char, &chr_index) == 4)
	while(1)
	{
		char* cur_line = getline(f_mrf);

		if(cur_line == NULL)
		{
			break;
		}

		if(sscanf(cur_line, "%s %c %d", mapping_map_str, &strand_char, &chr_index) != 3)
		{
			fprintf(stderr, "Could not parse fragment line: %s\n", cur_line);
		}

		int i_mapp_map = 0;
		t_string* cur_entry_length_str = new t_string();
		bool is_matching = false;
		char entry_type_char;

		//fprintf(stderr, "Processing cigar string: %s\n", quality_str);
		bool is_read_spliced = false;
		bool mapping_map_str_valid = validate_mapping_map_str(mapping_map_str, is_read_spliced);

		//int left_posn = chr_index;

		// If loading of all the fragments requested, process all the entries in the mapping string, then add all of them to the list of fragments for this read.
		while(mapping_map_str_valid && 
			mapping_map_str[i_mapp_map] != 0)
		{
			int l_cur_entry;
			get_next_entry_per_mapp_map_string(mapping_map_str,
												i_mapp_map, 
												is_matching,
												l_cur_entry,
												entry_type_char);

			// Analyze the fragment: Check the leading and following 'N's. This affects the length of the fragment.
			if(is_matching && per_fragment_callback != NULL)
			{
				// Call the fragment callback.
				per_fragment_callback(mapping_map_str, strand_char, chr_index, per_fragment_callback_param);
			}

			// Update the base for the current entry.
			// Must check whether to update the mapping posn: Update only for D and M entries.
			/*if(entry_type_char == 'D' || 
				entry_type_char == 'M' ||
				entry_type_char == 'N' ||
				entry_type_char == 'H')*/
			if(check_genome_index_update_per_CIGAR_entry(entry_type_char))
			{
				chr_index += l_cur_entry;
			}
		} // mapping map string processing loop.

		// Call the read callback.
		per_read_callback(mapping_map_str, strand_char, chr_index, per_read_callback_param);

		//fprintf(stderr, "%s: %s, %d (%d), %c\n", cur_line, new_mapped_read->mapping_str, new_mapped_read->base_index, new_mapped_read->span, new_mapped_read->strand);
		//getc(stdin);

		delete(cur_entry_length_str);
		delete [] cur_line;
	} // curent fragment data reading loop.

	close_f(f_mrf, mapped_reads_fp);
}

void add_mapped_fragments_per_mapped_read(t_mapped_read* mapped_read, vector<t_mapped_fragment*>* mapped_fragments)
{
	int i_mapp_map = 0;
	t_string* cur_entry_length_str = new t_string();
	bool is_matching = false;
	char entry_type_char;
	char strand_char = mapped_read->strand;
	//int chr_index = (mapped_read->strand=='F')?(mapped_read->base_index):(mapped_read->base_index-mapped_read->span+1);
	int chr_index = (mapped_read->base_index);
	char* mapping_map_str = mapped_read->mapping_str;

	//fprintf(stderr, "Processing cigar string: %s (%d, %c)\n", mapping_map_str, chr_index, strand_char);
	bool is_read_spliced = false;
	bool mapping_map_str_valid = validate_mapping_map_str(mapping_map_str, is_read_spliced);

	// If loading of all the fragments requested, process all the entries in the mapping string, then add all of them to the list of fragments for this read.
	while(mapping_map_str_valid && 
		mapping_map_str[i_mapp_map] != 0)
	{
		int l_cur_entry = 0;
		get_next_entry_per_mapp_map_string(mapping_map_str,
											i_mapp_map, 
											is_matching,
											l_cur_entry,
											entry_type_char);

		// Analyze the fragment: Check the leading and following 'N's. This affects the length of the fragment.		

		//int l_fragment = strlen(cur_fragment);
		if(is_matching)
		{
			//int l_fragment = get_l_fragment_per_cigar(quality_str);
			if(strand_char == 'F')
			{
				t_mapped_fragment* new_fragment = new t_mapped_fragment();
				new_fragment->base_index = chr_index;
				new_fragment->strand_char = strand_char;
				new_fragment->sequenced_fragment_length = l_cur_entry;
		
				mapped_fragments->push_back(new_fragment);

				//fprintf(stderr, "Adding: %c, %d (%d)\n", new_fragment->strand_char, new_fragment->base_index, new_fragment->sequenced_fragment_length);
				//getc(stdin);
			}
			else if(strand_char == 'R')
			{
				// Allocate and initialize a fragment and add it to the reverse strand fragment list.			
				t_mapped_fragment* new_fragment = new t_mapped_fragment();
				//new_fragment->base_index = chr_index + l_cur_entry - 1;
				new_fragment->base_index = chr_index;
				new_fragment->strand_char = strand_char;
				new_fragment->sequenced_fragment_length = l_cur_entry;

				mapped_fragments->push_back(new_fragment);

				//fprintf(stderr, "Adding: %c, %d (%d)\n", new_fragment->strand_char, new_fragment->base_index, new_fragment->sequenced_fragment_length);
				//getc(stdin);
				//rev_strand_frags->push_back(new_fragment);
			} // reverse strand check.
		} // maching check.

		// Update the base for the current entry.
		// Must check whether to update the mapping posn: Update only for D and M entries.
		//if(entry_type_char == 'D' || 
		//	entry_type_char == 'M' ||
		//	entry_type_char == 'N' ||
		//	entry_type_char == 'H')
		if(check_genome_index_update_per_CIGAR_entry(entry_type_char))
		{
			chr_index += l_cur_entry;
		}
	} // mapping map string processing loop.

	delete cur_entry_length_str;
}

// Following should be current with SAM specs.
bool check_read_nuc_index_update_per_CIGAR_entry(char entry_char)
{
	if(entry_char == '=' ||
		entry_char == 'X' ||
		entry_char == 'S' ||
		entry_char == 'M' ||
		entry_char == 'I')
	{
		return(true);
	}

	return(false);
}

// Following should be current with SAM specs.
bool check_genome_index_update_per_CIGAR_entry(char entry_char)
{
	if(entry_char == 'D' || 
		entry_char == 'M' ||
		entry_char == 'N' ||
		entry_char == '=' ||
		entry_char == 'X')
	{
		return(true);
	}

	return(false);
}

void exclude_reads_per_regions_per_chr(char* read_chr_id, 
	vector<t_mapped_read*>* cur_chr_reads,
	vector<t_mapped_read*>* no_overlap_reads,
	vector<t_annot_region*>* regions_2_exclude)
{
	if(regions_2_exclude == NULL)
	{
		for(int i_r = 0; i_r < (int)cur_chr_reads->size(); i_r++)
		{
			no_overlap_reads->push_back(cur_chr_reads->at(i_r));
		} // i_r loop.

		return;
	}

	// Restructure the regions.
	t_restr_annot_region_list* restructured_regions = restructure_annot_regions(regions_2_exclude);
	
	int n_total_excluded_reads = 0;
	
	int n_excluded_reads = 0;

	int i_reg_chr = t_string::get_i_str(restructured_regions->chr_ids, read_chr_id);

	if(i_reg_chr < (int)restructured_regions->chr_ids->size())
	{
		fprintf(stderr, "Excluding the reads in %s\n", restructured_regions->chr_ids->at(i_reg_chr));			

		for(int i_read = 0; i_read < (int)cur_chr_reads->size(); i_read++)
		{
			int cur_left_base_posn = (cur_chr_reads->at(i_read)->base_index);
			int cur_right_base_posn = (cur_chr_reads->at(i_read)->base_index + cur_chr_reads->at(i_read)->span - 1);

			int left_reg_i = locate_posn_per_sorted_obj_list(cur_left_base_posn, (vector<void*>*)(restructured_regions->regions_per_chrom[i_reg_chr]), 0, restructured_regions->regions_per_chrom[i_reg_chr]->size()-1, region_5p_accessor);

			bool cur_read_overlaps = false;
			while(left_reg_i > 0 &&
				restructured_regions->regions_per_chrom[i_reg_chr]->at(left_reg_i)->end > cur_left_base_posn)
			{
				left_reg_i--;
			}

			while(left_reg_i < (int)restructured_regions->regions_per_chrom[i_reg_chr]->size() &&
				restructured_regions->regions_per_chrom[i_reg_chr]->at(left_reg_i)->start < cur_right_base_posn)
			{
				int ol_start = MAX(restructured_regions->regions_per_chrom[i_reg_chr]->at(left_reg_i)->start, cur_left_base_posn);
				int ol_end = MIN(restructured_regions->regions_per_chrom[i_reg_chr]->at(left_reg_i)->end, cur_right_base_posn);

				if(ol_end >= ol_start)
				{
					cur_read_overlaps = true;
					break;
				}

				left_reg_i++;
			} // go over the region for the read.

			if(cur_read_overlaps)
			{
if(__DUMP_MAPPED_READ_TOOLS_MSGS__)
{
				if(cur_chr_reads->at(i_read)->strand == 'R')
				{
					fprintf(stderr, "Read: %s:%d(%d)(%c) overlaps with %s:%d-%d\n", read_chr_id, cur_chr_reads->at(i_read)->base_index, cur_chr_reads->at(i_read)->span, cur_chr_reads->at(i_read)->strand, 
						restructured_regions->chr_ids->at(i_reg_chr), restructured_regions->regions_per_chrom[i_reg_chr]->at(left_reg_i)->start, restructured_regions->regions_per_chrom[i_reg_chr]->at(left_reg_i)->end);
					getc(stdin);
				}
}

				// Update # of excluded reads.
				n_excluded_reads++;
			}
			else
			{
				no_overlap_reads->push_back(cur_chr_reads->at(i_read));
			}
		} // i_read loop.
	} // i_reg_chr search check.
	else
	{
		fprintf(stderr, "Skipping the regions in %s\n", read_chr_id);
		for(int i_read = 0; i_read < (int)cur_chr_reads->size(); i_read++)
		{
			no_overlap_reads->push_back(cur_chr_reads->at(i_read));
		} // i_read loop.
	} // read pile chromosome check.

	fprintf(stderr, "Excluded %d reads.\n", n_excluded_reads);

	n_total_excluded_reads += n_excluded_reads;

	delete_restructured_annot_regions(restructured_regions);

	fprintf(stderr, "Excluded %d reads in total.\n", n_total_excluded_reads);
}

void exclude_reads_per_regions(vector<char*>* read_chr_ids, 
	vector<vector<t_mapped_read*>*>* reads_per_chrs, 
	vector<vector<t_mapped_read*>*>* no_overlap_reads_per_chrs,
	vector<t_annot_region*>* regions_2_exclude)
{
	for(int i_read_chr = 0; i_read_chr < (int)read_chr_ids->size(); i_read_chr++)
	{
		vector<t_mapped_read*>* cur_chr_no_overlap_reads = new vector<t_mapped_read*>();
		no_overlap_reads_per_chrs->push_back(cur_chr_no_overlap_reads);

		exclude_reads_per_regions_per_chr(read_chr_ids->at(i_read_chr), 
			reads_per_chrs->at(i_read_chr),
			cur_chr_no_overlap_reads,
			regions_2_exclude);
	} // i_read_chr loop.
}

// Merge all the fragments of all the reads.
void get_mapped_fragments_per_mapped_reads(vector<t_mapped_read*>* mapped_reads, vector<t_mapped_fragment*>* mapped_fragments)
{
	for(int i_r = 0; i_r < (int)mapped_reads->size(); i_r++)
	{
		add_mapped_fragments_per_mapped_read(mapped_reads->at(i_r), mapped_fragments);
	} // i_r loop.
}

void delete_mapped_reads(vector<t_mapped_read*>* mapped_reads)
{
	for(int i_r = 0; i_r < (int)mapped_reads->size(); i_r++)
	{
		delete_mapped_read(mapped_reads->at(i_r));
	} // i_r loop.

	delete(mapped_reads);
}

void delete_mapped_read(t_mapped_read* mapped_read)
{
	delete [] mapped_read->mapping_str;
	delete(mapped_read);
}

/*
Prune reads:
Note that the trick here is to deal with the reads that have exact same pattern of mapping. We do not care about the
strand, this should be taken care of before the function is called.

Note that the pruning must be done at the read level, not at the fragment level.
*/
void prune_reads(vector<t_mapped_read*>* mapped_reads, int n_max_reps_per_posn, 
	vector<t_mapped_read*>* pruned_forward_reads, 
	vector<t_mapped_read*>* pruned_reverse_reads)
{
	// If the pruning is not requested, return all the reads.
	if(n_max_reps_per_posn == 0)
	{
		fprintf(stderr, "Skipping pruning.\n");
		for(int i_r = 0; i_r < (int)mapped_reads->size(); i_r++)
		{
			if(mapped_reads->at(i_r)->strand == 'F')
			{
				pruned_forward_reads->push_back(mapped_reads->at(i_r));
			}
			else
			{
				pruned_reverse_reads->push_back(mapped_reads->at(i_r));
			}
		} // i_r loop.

		return;
	}

	// Sort the mapped reads with respect to their 5' posn.
	sort(mapped_reads->begin(), mapped_reads->end(), sort_mapped_reads_per_5p);

    // First get rid of the extra fragments on forward strand.
	int* rep_cnts = new int[mapped_reads->size() + 2];
	memset(rep_cnts, 0, sizeof(int) * (mapped_reads->size() + 1));

	int prev_read_left_posn = 0;
    for(int i_r = 0; i_r < (int)mapped_reads->size(); i_r++)
    {
		int cur_read_left_posn = (mapped_reads->at(i_r)->base_index);
		
		if(i_r > 0 &&
			prev_read_left_posn == cur_read_left_posn)
		{
			rep_cnts[i_r] = rep_cnts[i_r-1] + 1;
		}
		else // This is a new fragment set its copy number to 1.
		{
			rep_cnts[i_r] = 1;
		}

		// Update prev. left posn.
		prev_read_left_posn = cur_read_left_posn;
    } // i_frag loop.

    for(int i_r = 0; i_r < (int)mapped_reads->size(); i_r++)
    {
		if(rep_cnts[i_r] <= n_max_reps_per_posn)
		{
			if(mapped_reads->at(i_r)->strand == 'F')
			{
				pruned_forward_reads->push_back(mapped_reads->at(i_r));
			}
			else
			{
				pruned_reverse_reads->push_back(mapped_reads->at(i_r));
			}
		}
		else // This is a new fragment set its copy number to 1.
		{
			// Delete the pruned reads, otherwise they will be lost.
			delete_mapped_read(mapped_reads->at(i_r));

if(__DUMP_MAPPED_READ_TOOLS_MSGS__)
{
			fprintf(stderr, "Pruning %d repetition read @ %d, %c.\n", rep_cnts[i_r], mapped_reads->at(i_r)->base_index, mapped_reads->at(i_r)->strand);
			getc(stdin);
}
		}
    } // i_r loop.

	delete [] rep_cnts;

if(__DUMP_MAPPED_READ_TOOLS_MSGS__)
	fprintf(stderr, "Pruned to %ld forward, %ld reverse strand fragments.\n", pruned_forward_reads->size(), pruned_reverse_reads->size());
}

bool sort_mapped_reads_per_5p(t_mapped_read* read1, t_mapped_read* read2)
{
	int frag1_5p = read_5p_accessor(read1);
	int frag2_5p = read_5p_accessor(read2);

	return(frag1_5p < frag2_5p);
}

int read_5p_accessor(void* obj_ptr)
{
	t_mapped_read* frag_obj_ptr = (t_mapped_read*)obj_ptr;

	return(frag_obj_ptr->base_index);
}

//int read_3p_accessor(void* obj_ptr)
//{
//	t_mapped_read* frag_obj_ptr = (t_mapped_read*)obj_ptr;
//
//	if(frag_obj_ptr->strand == 'F')
//	{
//		return(frag_obj_ptr->base_index+frag_obj_ptr->span-1);
//	}
//	else if(frag_obj_ptr->strand == 'R')
//	{
//		return(frag_obj_ptr->base_index);
//	}
//	else
//	{
//		fprintf(stderr, "The strand char for fragment object is %s @ %s(%d)\n", frag_obj_ptr->strand, __FILE__, __LINE__);
//		exit(0);
//		return(-1);
//	}
//}

