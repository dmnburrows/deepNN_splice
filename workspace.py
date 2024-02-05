#!/bin/bash
# filter reads quality? // trim_galore
%%bash
path=/cndd3/dburrows/DATA/splice/smartseq-rna_MOp_biccn/pseudobulk/coarse/
#id_arr=($(find $inpath/ -maxdepth 1 -name "*Sample*$str*"))
cell_arr=("GLU" "GABA")
for c in ${cell_arr[@]}
do
    curr=($(ls $path$c))
    echo ${curr[0]}
    echo ${curr[1]}
    trim_galore --cores 10 --paired -q 20 -e 0.1 --fastqc_args "-noextract" $path$c/${curr[0]} $path$c/${curr[1]} -o $path$c
    
done
echo Done