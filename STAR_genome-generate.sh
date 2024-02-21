#!/bin/bash
echo Starting
STAR --runMode genomeGenerate \
     --runThreadN 12 \
     --genomeDir /cndd3/dburrows/DATA/annotations/genome/grcm38.p3/STAR_genome.1-19-X \
     --genomeFastaFiles /cndd3/dburrows/DATA/annotations/genome/grcm38.p3/GRCm38.p3.genome_1-19-X.fa \
     --sjdbGTFfile /cndd3/dburrows/DATA/annotations/genome/grcm38.p3/gencode.vM3.chr_1-19-X.annotation.gtf \
     --sjdbOverhang 50
echo DoneS