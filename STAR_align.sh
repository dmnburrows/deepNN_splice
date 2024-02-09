#!/bin/bash
#Define variables
r1=$1 #r1 path
r2=$2 #r2 path
out=$3 #out path



STAR --genomeDir /cndd3/dburrows/DATA/annotations/genome/grcm38.p3/STAR_genome \
     --readFilesIn $r1 $r2 \
     --runThreadN 14 \
     --outFileNamePrefix $out \
     --readFilesCommand zcat \
     --sjdbGTFfile /cndd3/dburrows/DATA/annotations/genome/grcm38.p3/gencode.vM3.chr_patch_hapl_scaff.annotation.gtf \
     --outSAMtype BAM SortedByCoordinate \
     --sjdbOverhang 50 \
     --alignSJoverhangMin 8 \
     --alignSJDBoverhangMin 3 \
     --alignIntronMin 20 \
     --alignIntronMax 1000000 \
     --alignMatesGapMax 1000000 \
     --outFilterMismatchNmax 999 \
     --outFilterMultimapNmax 20 \
     --outFilterType BySJout \
     --outSAMattributes All \
     --sjdbScore 1 \
     --genomeLoad LoadAndKeep \
     --limitBAMsortRAM 10000000000 \
     --quantMode TranscriptomeSAM GeneCounts \
     --winAnchorMultimapNmax 200 \
     --outMultimapperOrder Random \
     --outSAMmultNmax -1 \
     
     