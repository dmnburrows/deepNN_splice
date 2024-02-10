#!/bin/bash
#Run allcools merge to sum over regions
path=$1
out=$2
allcools merge-allc \
    --cpu 12 \
    --allc_paths $path \
    --output_path $out \
    --chrom_size_path /cndd3/dburrows/DATA/annotations/genome/grcm38.p3/modified_mm10.chrom.sizes
