#!/bin/bash
PATH_TO_OUTS="/home/vsevim/prj/1012-ckd/S1/analysis/primary/cellranger/Screen1_66guides_L1/outs"
PATH_TO_BC="${PATH_TO_OUTS}/filtered_feature_bc_matrix"
PATH_TO_REF="/home/vsevim/prj/refs/refdata-gex-GRCh38-2020-A/fasta/genome.fa"
num_threads_to_use=40
output_dir_name=soupor_L1
num_clusters=4

singularity exec \
	~/software/souporcell_latest.sif \
	souporcell_pipeline.py -i "${PATH_TO_OUTS}/possorted_genome_bam.bam" \
	-b "${PATH_TO_BC}/barcodes.tsv.gz" \
	-f "${PATH_TO_REF}" \
	-t ${num_threads_to_use} \
	-o ${output_dir_name} \
	-k ${num_clusters}