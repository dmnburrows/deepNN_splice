#!/bin/bash
#Define variables
r1=$1 #r1 path
r2=$2 #r2 path
out=$3 #out path



STAR --genomeDir /cndd3/dburrows/DATA/annotations/genome/grcm38.p3/STAR_genome.1-19-X \
     --readFilesIn $r1 $r2 \
     --runThreadN 15 \
     --outFileNamePrefix $out \
     --readFilesCommand zcat \
     --sjdbGTFfile /cndd3/dburrows/DATA/annotations/genome/grcm38.p3/gencode.vM3.chr_1-19-X.annotation.gtf \
     --outSAMtype BAM Unsorted \
     --sjdbOverhang 50 \
     --alignSJoverhangMin 20 \
     --alignSJDBoverhangMin 3 \
     --alignIntronMin 20 \
     --alignIntronMax 1000000 \
     --alignMatesGapMax 1000000 \
     --outFilterMismatchNmax 999 \
     --outFilterMultimapNmax 20 \
     --outFilterType BySJout \
     --outSAMattributes All \
     --sjdbScore 1 \
     --genomeLoad NoSharedMemory \
     --limitBAMsortRAM 85000000000 \
     --quantMode TranscriptomeSAM GeneCounts \
     --winAnchorMultimapNmax 200 \
     --outMultimapperOrder Random \
     --outSAMmultNmax -1 \
     --limitOutSJcollapsed 2000000 \
     --outSJfilterOverhangMin -1 -1 -1 -1 
     


     
