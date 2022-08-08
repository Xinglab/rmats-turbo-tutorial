rmats2sashimiplot  \
	--b1 ../bam/PC3E_rep1/Aligned.sortedByCoord.out.bam,../bam/PC3E_rep2/Aligned.sortedByCoord.out.bam,../bam/PC3E_rep3/Aligned.sortedByCoord.out.bam \
	--b2 ../bam/GS689_rep1/Aligned.sortedByCoord.out.bam,../bam/GS689_rep2/Aligned.sortedByCoord.out.bam,../bam/GS689_rep3/Aligned.sortedByCoord.out.bam \
	-t SE \
	-e sashimi_events_example.txt \
	--l1 PC3E_rep --l2 GS689_rep \
	--exon_s 1 --intron_s 5 \
	-o ./output_example \
	--group-info sashimi_groupInfo.txt
