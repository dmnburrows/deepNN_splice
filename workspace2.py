#!/bin/bash
# sort + index
path=/cndd2/dburrows/DATA/splice/smartseq-rna_MOp_biccn/pseudobulk/gran1/
code_path=/cndd3/dburrows/CODE/deepNN_splice/workspace2.py

cell_arr=('Vip' 'Pvalb' 'Sst' 'Lamp5' )
for c in ${cell_arr[@]}
do
    echo $path$c/Aligned.out.bam
    samtools sort -o $path$c/Aligned.sortedByCoord.out.bam -m 5G -@ 10 $path$c/Aligned.out.bam
    samtools index $path$c/Aligned.sortedByCoord.out.bam
    
done
echo Done

