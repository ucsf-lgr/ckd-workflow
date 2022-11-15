#!/bin/bash
set -euo pipefail

txtred=$'\e[0;31m' # Red
txtgrn=$'\e[0;32m' # Green
txtylw=$'\e[0;33m' # Yellow
txtwht=$'\e[0;37m'

# prj_name  <- "CKD_Screen1_Lib1_66guides"
# prj_path  <- "/home/vsevim/prj/1012-ckd/S1/"
# data_subdir <- "analysis/primary/cellranger/Screen1_66guides_L1/outs"

SAVE_H5="YES"

for LIB_NO in 1 2 3 4; do
    #LIB_NO="1"
    PRJ_NAME="Screen1_66guides"
    LIBRARY_NAME="Lib_${LIB_NO}"
    H5_FILE_NAME="${PRJ_NAME}_${LIBRARY_NAME}.h5seurat"

    PRJ_PATH="/home/vsevim/prj/1012-ckd/S1/"
    CELLRANGER_OUT_SUBDIR="/analysis/primary/cellranger/${PRJ_NAME}_L${LIB_NO}/outs/"

    SECONDARY_A_SCRIPT_PATH="/home/vsevim/prj/workflows/ckd/secondary/"
    PRIMARY_A_OUT_PATH="${PRJ_PATH}/analysis/primary/${LIBRARY_NAME}/"
    SECONDARY_A_OUT_PATH="${PRJ_PATH}/analysis/secondary/${LIBRARY_NAME}/"
    SEURAT_OBJ_PATH="$SECONDARY_A_OUT_PATH/seurat_objects/"
    NOTEBOOK_OUT_PATH="$SECONDARY_A_OUT_PATH/notebooks/"

    mkdir -p ${NOTEBOOK_OUT_PATH}
    mkdir -p ${SEURAT_OBJ_PATH}

    # STEP 1 ----------------------------------------------------------------------------
    IPYNB_FILE="00_QC.R.ipynb"
    OUT_IPYNB="00_QC__${PRJ_NAME}_${LIBRARY_NAME}.ipynb"
    JUPYTER_FILE_TO_EXECUTE="$SECONDARY_A_SCRIPT_PATH/$IPYNB_FILE"
    JUPYTER_FILE_TO_PRODUCE="${NOTEBOOK_OUT_PATH}/${OUT_IPYNB}"

    echo -e "$txtgrn \n*** STEP 1: $IPYNB_FILE $txtwht"
    echo -e "$txtgrn ${JUPYTER_FILE_TO_EXECUTE} $txtwht"
    echo -e ${SECONDARY_A_SCRIPT_PATH}
    echo -e ${PRIMARY_A_OUT_PATH}
    echo -e ${SECONDARY_A_OUT_PATH}
    echo -e ${NOTEBOOK_OUT_PATH}
    echo -e ${SEURAT_OBJ_PATH}/${H5_FILE_NAME}


    papermill \
        -p prj_name ${PRJ_NAME}    \
        -p prj_path ${PRJ_PATH}    \
        -p data_subdir ${CELLRANGER_OUT_SUBDIR}  \
        -p seurat_obj_path ${SEURAT_OBJ_PATH} \
        -p seurat_obj_fname ${H5_FILE_NAME} \
        -p library_name ${LIBRARY_NAME} \
        -p save_seurat_h5 ${SAVE_H5} \
        -p papermill_run "YES" \
        "${JUPYTER_FILE_TO_EXECUTE}" \
        "${JUPYTER_FILE_TO_PRODUCE}"

    jupyter nbconvert --to html --no-input --no-prompt ${JUPYTER_FILE_TO_PRODUCE}
done