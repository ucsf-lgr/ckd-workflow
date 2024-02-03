#!/bin/bash

#echo "Update the dirs etc to run this script!"
#exit -1

TRANSCRIPTOME="refs/refdata-gex-GRCh38-2020-A"
FEATURE_REF="S1_resources/guide_reference.csv"
cellranger=cellranger/cellranger-7.0.0/cellranger

cellranger_out_path="S1/analysis/cellranger"
log_path="S1/scripts/logs"
config_path="S1/scripts/library_csv"
decsription="CKD project Screen 1 perturb-seq run with 
             66 guides, 4 of which are NT controls."


SCREEN="Screen1_66guides"
LOG="$ID.log"
ERR="$ID.err"
        
for ID in L1 L2 L3 L4
do
        LOG="$ID.log"
        ERR="$ID.err"
	SCR_NAME="${SCREEN}_${ID}"
        LIBRARY_REF="$config_path/${ID}_library_reference.csv"
        cd "$cellranger_out_path"
        echo "Running $LIBRARY_REF"

        "$cellranger" count \
                --id="$SCR_NAME" \
                --description="$decription" \
                --libraries="$LIBRARY_REF" \
                --transcriptome="$TRANSCRIPTOME" \
                --feature-ref="$config_path/$FEATURE_REF" \
                --localcores=40 \
                --localmem=400 \
                --no-bam \
                --nosecondary \
                --disable-ui \
                --expect-cells=10000 \
                > "$log_path/$LOG" 2> "$log_path/$ERR"                
done >"$log_path/00.stdout" 2>"$log_path/00.stdout" 
