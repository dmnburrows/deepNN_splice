#!/bin/bash
# filter reads quality? // trim_galore
path=/cndd3/dburrows/DATA/splice/smartseq-rna_MOp_biccn/pseudobulk/gran1/
scp /cndd3/dburrows/CODE/deepNN_splice/workspace.py $path/log.workspace.trim_galore
cell_arr=("Vip" "Pvalb" "Sst" "Lamp5")

for c in ${cell_arr[@]}
do
    curr=($(ls $path$c/*fastq.gz ))
    echo ${curr[0]}
    echo ${curr[1]}
    trim_galore --cores 10 --paired -q 20 -e 0.1 --fastqc_args "-noextract" ${curr[0]} ${curr[1]} -o $path$c
    
done
echo Done
