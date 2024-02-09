#!/bin/bash
echo Starting
STAR --runMode genomeGenerate \
     --runThreadN 12 \
     --genomeDir /cndd3/dburrows/DATA/annotations/genome/grcm38.p3 \
     --genomeFastaFiles /cndd3/dburrows/DATA/annotations/genome/grcm38.p3/GRCm38.p3.genome.fa.gz \
     --sjdbGTFfile /cndd3/dburrows/DATA/annotations/genome/grcm38.p3/gencode.vM3.chr_patch_hapl_scaff.annotation.gtf.gz \
     --sjdbOverhang 50
echo Done