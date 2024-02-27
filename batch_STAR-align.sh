#!/bin/bash
# align with STAR
in_path=/cndd3/dburrows/DATA/splice/smartseq-rna_MOp_biccn/pseudobulk/gran1/
out_path=/cndd2/dburrows/DATA/splice/smartseq-rna_MOp_biccn/pseudobulk/gran1/
code_path=/cndd3/dburrows/CODE/deepNN_splice/STAR_align.sh
chmod u+x $code_path
scp $code_path $out_path/log.workspace.STAR_align

cell_arr=("L6_CT")
for c in ${cell_arr[@]}
do
    curr=($(ls $in_path$c/*val*fq.gz))
    echo ${curr[0]}
    echo ${curr[1]}
    $code_path ${curr[0]} ${curr[1]} $out_path$c/
    #samtools index $out_path$c/Aligned.sortedByCoord.out.bam
    
done
echo Done