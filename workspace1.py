#!/bin/bash
# filter reads quality? // trim_galore
path=/cndd3/dburrows/DATA/splice/smartseq-rna_MOp_biccn/pseudobulk/gran1/
scp /cndd3/dburrows/CODE/deepNN_splice/workspace.py $path/log.workspace.trim_galore
cell_arr=("L2-3_IT" "L5_IT" "L6b" "L6_CT" "L6_IT")

for c in ${cell_arr[@]}
do
    curr=($(ls $path$c/*fastq.gz ))
    echo ${curr[0]}
    echo ${curr[1]}
    trim_galore --cores 10 --paired -q 20 -e 0.1 --fastqc_args "-noextract" ${curr[0]} ${curr[1]} -o $path$c
    
done
echo Done


#!/bin/bash
# align with STAR

# path=/cndd3/dburrows/DATA/splice/smartseq-rna_MOp_biccn/pseudobulk/coarse/
# code_path=/cndd3/dburrows/CODE/deepNN_splice/STAR_align.sh
# chmod u+x $code_path
# scp $code_path $path/log.workspace.STAR_align

# cell_arr=("GLU" "GABA")
# for c in ${cell_arr[@]}
# do
#     curr=($(ls $path$c/*val*fq.gz))
#     echo ${curr[0]}
#     echo ${curr[1]}
#     $code_path ${curr[0]} ${curr[1]} $path$c/
#     samtools index $path$c/Aligned.sortedByCoord.out.bam
    
# done
# echo Done