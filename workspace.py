#!/bin/bash
# align with STAR

path=/cndd3/dburrows/DATA/splice/smartseq-rna_MOp_biccn/pseudobulk/coarse/
code_path=/cndd3/dburrows/CODE/deepNN_splice/STAR_align.sh
chmod u+x $code_path
scp $code_path $path/log.workspace

cell_arr=("GLU" "GABA")
for c in ${cell_arr[@]}
do
    curr=($(ls $path$c/*val*fq.gz))
    echo ${curr[0]}
    echo ${curr[1]}
    ./$code_path ${curr[0]} ${curr[1]} $path$c/
    samtools index $path$c/Aligned.sortedByCoord.out.bam
    
done
echo Done