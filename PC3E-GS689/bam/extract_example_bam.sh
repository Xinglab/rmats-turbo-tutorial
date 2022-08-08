samtools view -b /mnt/isilon/xing_lab/wangy14/project/annotation_rmats-turbo/PC3E-GS689/bam/GS689_rep3/Aligned.sortedByCoord.out.bam chr4:76649777-76735439 > USO1_Aligned.sortedByCoord.out.bam 

samtools view -b /mnt/isilon/xing_lab/wangy14/project/annotation_rmats-turbo/PC3E-GS689/bam/GS689_rep3/Aligned.sortedByCoord.out.bam chr19:18208603-18262502 > MAST3_Aligned.sortedByCoord.out.bam

rm test.bam 2> /dev/null
samtools merge test.bam MAST3_Aligned.sortedByCoord.out.bam USO1_Aligned.sortedByCoord.out.bam
samtools sort test.bam -o Aligned.sortedByCoord.out.bam
rm test.bam
rm USO1_Aligned.sortedByCoord.out.bam
rm MAST3_Aligned.sortedByCoord.out.bam

