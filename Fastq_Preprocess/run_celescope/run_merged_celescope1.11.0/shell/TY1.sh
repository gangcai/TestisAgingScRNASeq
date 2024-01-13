celescope rna sample --outdir TY1/TY1/00.sample --sample TY1 --thread 15 --chemistry auto  --fq1 ../transcriptome_seq/P22033104/rawdata/TY1_R1.fq.gz 
celescope rna barcode --outdir TY1/TY1/01.barcode --sample TY1 --thread 15 --chemistry auto --lowNum 2 --allowNoPolyT  --fq1 ../transcriptome_seq/P22033104/rawdata/TY1_R1.fq.gz --fq2 ../transcriptome_seq/P22033104/rawdata/TY1_R2.fq.gz 
celescope rna cutadapt --outdir TY1/TY1/02.cutadapt --sample TY1 --thread 15 --minimum_length 20 --nextseq_trim 20 --overlap 10 --insert 150  --fq TY1/TY1/01.barcode/TY1_2.fq 
celescope rna star --outdir TY1/TY1/03.star --sample TY1 --thread 15 --genomeDir /home/db/public/SoftwareIndex/celescope1.11.0/mouse --outFilterMultimapNmax 1 --starMem 30  --fq TY1/TY1/02.cutadapt/TY1_clean_2.fq 
celescope rna featureCounts --outdir TY1/TY1/04.featureCounts --sample TY1 --thread 15 --gtf_type exon --genomeDir /home/db/public/SoftwareIndex/celescope1.11.0/mouse  --input TY1/TY1/03.star/TY1_Aligned.sortedByCoord.out.bam 
celescope rna count --outdir TY1/TY1/05.count --sample TY1 --thread 15 --genomeDir /home/db/public/SoftwareIndex/celescope1.11.0/mouse --expected_cell_num 50000 --cell_calling_method EmptyDrops_CR  --bam TY1/TY1/04.featureCounts/TY1_name_sorted.bam --force_cell_num None 
celescope rna analysis --outdir TY1/TY1/06.analysis --sample TY1 --thread 15 --genomeDir /home/db/public/SoftwareIndex/celescope1.11.0/mouse  --matrix_file TY1/TY1/05.count/TY1_filtered_feature_bc_matrix 
