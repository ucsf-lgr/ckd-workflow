#!/bin/bash

#echo "Update the dirs etc to run this script!"
#exit -1

TRANSCRIPTOME="/home/vsevim/prj/refs/refdata-gex-GRCh38-2020-A"
FEATURE_REF="guide_reference.csv"
cellranger=/home/sfederman/programs/cellranger/cellranger-7.0.0/cellranger

cellranger_out_path="/home/vsevim/prj/1012-ckd/S1/analysis/cellranger"
data_path="/home/vsevim/prj/1012-ckd/S1/data/"
script_path="/home/vsevim/prj/1012-ckd/S1/scripts"
log_path="/home/vsevim/prj/1012-ckd/S1/scripts/logs"
config_path="/home/vsevim/prj/1012-ckd/S1/scripts/library_csv"
decsription="CKD project Screen 1perturb-seq run with 66 guides, 
             4 of which are NT controls. Nov 7 2022"


SCREEN="Screen1_66guides"
LOG="$ID.log"
ERR="$ID.err"

for ID in L1 #L2 L3 L4
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
                --disable-ui \
                --expect-cells=10000 \
                > "$log_path/$LOG" 2> "$log_path/$ERR"                
done >"$log_path/00.stdout" 2>"$log_path/00.stdout" 
