FILE_INIT=1
ANNOTATION_GFF_GZ_FP=gencode.v31.annotation.gff3.gz
ASSEMBLYID=hg38
GENOME_DIR=/internal/aharmanci1/dir/genomes/hg38
SAMTOOLS_PATH=samtools
PC_GENE_NAME_ID_MAP=PC_gene_id_name.map

if [ ${FILE_INIT} == 1 ]
then
        wget -c ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_31/gencode.v31.annotation.gff3.gz

        # Extract the mapper.
        gzip -cd ${ANNOTATION_GFF_GZ_FP} | grep "gene_type=protein_coding;" | awk {'if($3=="gene"){info=$9;split(info, arr, ";");for(i=1;i<length(arr);i++){split(arr[i], arr2, "=");if(arr2[1]=="gene_name"){gene_name=arr2[2]};if(arr2[1]=="gene_id"){gene_id=arr2[2]}};split(gene_id, arr3, ".");print gene_name"\t"arr3[1]}'} | sort -u -k1,1 | sort -u -k2,2 > ${PC_GENE_NAME_ID_MAP}

        # Write the tx2gene id mapping ile.
        gzip -cd ${ANNOTATION_GFF_GZ_FP} | awk 'BEGIN{print "TXNAME\tGENEID"}{if($3=="transcript"){info=$9;split(info, arr, ";");for(i=1;i<length(arr);i++){split(arr[i], arr2, "=");if(arr2[1]=="transcript_id"){transcript_id=arr2[2]};if(arr2[1]=="gene_id"){gene_id=arr2[2]}};if(length(transcript_id)>0 && length(gene_id)>0){print transcript_id"\t"gene_id}}}' > tx2gene.txt

        fetchChromSizes ${ASSEMBLYID} > ${ASSEMBLYID}.list
        sed -i 's/chr//g' ${ASSEMBLYID}.list
fi

RUN_GET_GENES_INTERVAL=0
if [ $RUN_GET_GENES_INTERVAL == 1 ]
then
        # Extract the exons.
        gzip -cd ${ANNOTATION_GFF_GZ_FP} | grep -v "tag=PAR" | grep "gene_type=protein_coding;" | awk {'if($1!="chrY" && $1!="chrM" && $3=="exon")print $0'} > annotation_exons.gff
        ./bin/IntrExtract -GENCODE_GTF_2_Interval_per_feature annotation_exons.gff exon gene_id regions.interval
		
		# ABOVE COMMAND CREATES THE FILE CALLED regions.interval.
fi

BAM_fp=/internal/aharmanci1/dir/Klisch_Patel_MeningioMine_7.24.2019/SingleCell/2Samples/Front_analysis/outs/possorted_genome_bam.bam
barcodes_fp=/internal/aharmanci1/dir/Klisch_Patel_MeningioMine_7.24.2019/SingleCell/2Samples/Front_analysis/outs/filtered_feature_bc_matrix/barcodes.tsv.gz

./bin/IntrExtract -compute_single_cell_expression_stats_per_10X_SAM ${barcodes_fp} ${BAM_fp} regions.interval intrex_stats_ST.txt.gz
