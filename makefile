all: IntrExtract 

CC = g++
comp_flags = -c -Wall -O3
exec_name = bin/IntrExtract
LIB_DIR = src

# Define pattern rule for building object files.
%.o: %.cpp
	${CC} ${comp_flags} $< -o $@

objs = \
${LIB_DIR}/main.o \
${LIB_DIR}/ntrx_config.o \
${LIB_DIR}/ntrx_ansi_cli.o \
${LIB_DIR}/ntrx_nucleotide.o \
${LIB_DIR}/ntrx_ansi_string.o \
${LIB_DIR}/ntrx_genomics_coords.o \
${LIB_DIR}/ntrx_mapped_read_tools.o \
${LIB_DIR}/ntrx_genome_sequence_tools.o \
${LIB_DIR}/ntrx_signal_track_tools.o \
${LIB_DIR}/ntrx_expression_tools.o \
${LIB_DIR}/ntrx_nomenclature.o \
${LIB_DIR}/ntrx_annot_region_tools.o \
${LIB_DIR}/ntrx_gff_utils.o \
${LIB_DIR}/ntrx_human_data_processing.o \
${LIB_DIR}/ntrx_file_utils.o \
${LIB_DIR}/ntrx_histogram.o \

IntrExtract: ${objs}
	${CC} -O3 -lz -o ${exec_name} ${objs}

clean:
	rm -f ${objs} ${exec_name} 
