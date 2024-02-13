#!/bin/bash
echo Starting
STAR --runMode genomeGenerate \
     --runThreadN 15 \
     --genomeDir /cndd3/dburrows/DATA/annotations/genome/grcm38.p3/STAR_genome \
     --genomeFastaFiles /cndd3/dburrows/DATA/annotations/genome/grcm38.p3/GRCm38.p3.genome.fa \
     --sjdbGTFfile /cndd3/dburrows/DATA/annotations/genome/grcm38.p3/gencode.vM3.chr_patch_hapl_scaff.annotation.gtf \
     --sjdbOverhang 50
echo DoneS