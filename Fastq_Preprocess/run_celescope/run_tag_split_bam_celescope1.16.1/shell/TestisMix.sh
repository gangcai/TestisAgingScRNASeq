r1file="../transcriptome_seq/P22033104/rawdata/TY1_R1.fq.gz"
bam_file="/home/db/private/XieLab/aging/mouse_testis_scRNASeq/singleron_6samples/run_merged_celescope1.11.0/TestisMix/TestisMix/04.featureCounts/TestisMix_Aligned.sortedByCoord.out.bam.featureCounts.bam"
celescope tag split_tag --outdir .//TestisMix/06.split_tag --sample TestisMix --thread 15 --split_bam  --match_dir ../run_merged_celescope1.11.0/TestisMix/TestisMix/ --umi_tag_file .//TestisMix/04.count_tag/TestisMix_umi_tag.tsv --R1_read ${r1file} --bam_file ${bam_file}
