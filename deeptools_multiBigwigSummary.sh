#!/bin/bash

#### Merge bw files with deeptools ####
path=/cndd2/dburrows/DATA/splice/snmcseq-mc_MOp_biccn/pseudobulk/gran1/
data_l=("ATAC")
cell_l=($(ls $path))

for c in ${cell_l[@]}
do
    echo Running $c
    
    for d in ${data_l[@]}
    do
    echo Running $d
        bw_l=($(ls $path$c/*$d*bw))
        echo ${bw_l[@]}
        
        # Running multiBigwigSummary
        multiBigwigSummary bins \
            --bwfiles ${bw_l[@]} \
            --binSize 128 \
            --numberOfProcessors 12 \
            --outFileName $path$c/$c'_'$d'.merge.bingraph.npz' \
            --outRawCounts $path$c/$c'_'$d'.merge.bingraph.raw'
            echo "multiBigwigSummary analysis is complete."
    done
done
echo 'All Done'

