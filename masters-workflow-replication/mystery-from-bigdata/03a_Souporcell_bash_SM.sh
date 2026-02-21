#!/bin/bash

set -euo pipefail

fa=/lustre/scratch117/cellgen/cellgeni/TIC-starsolo/tic-1262/GRCh38_v32_modified.fa
im=/nfs/cellgeni/singularity/images/souporcell.sif
cv=/lustre/scratch117/cellgen/cellgeni/simon/altered_vcfs/chr_filtered_2p_1kgenomes_GRCh38.vcf

CPU=${1?Need CPU count}
dir=${2?Need directory}
thek=${3?Need k}

cd ../data/$dir
psb=removed.bam
bcd=barcodes.tsv

[[ ! -e "$psb" ]] && echo "No $psb in $dir" && false
[[ ! -e "$bcd" ]] && echo "No $bcd in $dir" && false

echo $dir ok
mkdir -p soc

singularity exec -B /lustre -B /nfs \
      $im souporcell_pipeline.py \
      -i $psb                 \
      -b $bcd                 \
      -f $fa                  \
      -k $thek                \
      --common_variants $cv   \
      -t $CPU -o soc        \
      --skip_remap True     \
      --no_umi True
