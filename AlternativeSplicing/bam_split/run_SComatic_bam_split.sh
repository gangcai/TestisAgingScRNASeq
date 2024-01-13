#!/bin/bash
root="/public/packages/SComatic/SComatic/"
script=${root}/scripts/SplitBam/SplitBamCellTypes.py
bam_dir="/home/db/private/XieLab/aging/mouse_testis_scRNASeq/singleron_6samples/clean/bam"
#sample="TO1"
for sample in TO1 TO2 TO3 TY1 TY2 TY3
do
	mkdir -p ${sample}
	bam=${bam_dir}/${sample}.bam
	meta=${sample}_cell_cluster_info.tsv
	python ${script} --bam $bam --meta $meta --id $sample --outdir $sample
done
