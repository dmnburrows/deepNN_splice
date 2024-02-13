#!/bin/bash
# align with STAR

outpath=/cndd3/dburrows/DATA/splice/snmcseq-mc_MOp_biccn/allc/pseudobulk/coarse/
inpath=/cndd3/dburrows/DATA/splice/snmcseq-mc_MOp_biccn/allc/
code_path=/cndd3/dburrows/CODE/deepNN_splice/workspace.py
chmod u+x $code_path
scp $code_path $outpath/log.workspace

cell_arr=("GLU" "GABA")
for c in "${cell_arr[@]}"
do
    echo $c
    for num in 1 2
    do
        echo Doing $num
        allcools merge-allc \
            --cpu 30 \
            --allc_paths $inpath$c/list$num.txt \
            --output_path $outpath/$c/$num/ \
            --chrom_size_path /cndd3/dburrows/DATA/annotations/genome/grcm38.p3/modified_mm10.chrom.sizes
        echo Done $num
    done
done
echo Done


