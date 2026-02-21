#! /bin/bash

set -euo pipefail

TAG=$1
WDIR=`pwd`
EPOCHS=150
EXPECTED_CELLS=10000
DROPLETS=20000
FPR=0.01
LEARN=0.0001

# add singularity to the path
PATH="/software/singularity-v3.6.4/bin:$PATH"

# set the path to the image we want to use
CELLBENDER_IMAGE="/nfs/cellgeni/singularity/images/cellbender-0.2.0.sif"

# path to output folder (samples will have a folder inside this one)
OUTPUT_FOLDER="$WDIR/results/$TAG/outs"

# create output folder if it does not exist
if [[ ! -d "${OUTPUT_FOLDER}" ]]; then
    mkdir -p "${OUTPUT_FOLDER}"
fi

INPUT_FILE="$WDIR/data/$TAG/raw_gene_bc_matrices_h5.h5"

singularity run --nv --bind /nfs,/lustre $CELLBENDER_IMAGE cellbender remove-background \
    --input "$INPUT_FILE" \
    --output "${OUTPUT_FOLDER}/cellbender_out.h5" \
    --cuda \
    --expected-cells $EXPECTED_CELLS \
    --epochs $EPOCHS \
    --total-droplets-included $DROPLETS \
    --fpr $FPR \
    --learning-rate $LEARN
