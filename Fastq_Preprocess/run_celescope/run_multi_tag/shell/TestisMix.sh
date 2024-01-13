celescope tag sample --outdir .//TestisMix/00.sample --sample TestisMix --thread 15 --chemistry auto  --fq1 ../tag_seq/P22033104/rawdata/TY1_R1.fq.gz 
celescope tag barcode --outdir .//TestisMix/01.barcode --sample TestisMix --thread 15 --chemistry auto --lowNum 2  --fq1 ../tag_seq/P22033104/rawdata/TY1_R1.fq.gz --fq2 ../tag_seq/P22033104/rawdata/TY1_R2.fq.gz 
celescope tag cutadapt --outdir .//TestisMix/02.cutadapt --sample TestisMix --thread 15 --minimum_length 20 --nextseq_trim 20 --overlap 10 --insert 150  --fq .//TestisMix/01.barcode/TestisMix_2.fq 
celescope tag mapping_tag --outdir .//TestisMix/03.mapping_tag --sample TestisMix --thread 15 --fq_pattern L25C15 --barcode_fasta ../Clindex/tag_barcode.fasta  --fq .//TestisMix/02.cutadapt/TestisMix_clean_2.fq 
celescope tag count_tag --outdir .//TestisMix/04.count_tag --sample TestisMix --thread 15 --UMI_min auto --dim 1 --SNR_min auto --coefficient 0.1  --match_dir ../run_merged_celescope1.11.0/TestisMix/TestisMix/ --read_count_file .//TestisMix/03.mapping_tag/TestisMix_read_count.tsv 
celescope tag analysis_tag --outdir .//TestisMix/05.analysis_tag --sample TestisMix --thread 15  --match_dir ../run_merged_celescope1.11.0/TestisMix/TestisMix/ --tsne_tag_file .//TestisMix/04.count_tag/TestisMix_tsne_tag.tsv 
celescope tag split_tag --outdir .//TestisMix/06.split_tag --sample TestisMix --thread 15 --split_fastq --split_matrix  --match_dir ../run_merged_celescope1.11.0/TestisMix/TestisMix/ --umi_tag_file .//TestisMix/04.count_tag/TestisMix_umi_tag.tsv 
