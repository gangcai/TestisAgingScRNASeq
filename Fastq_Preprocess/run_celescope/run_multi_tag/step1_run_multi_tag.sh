#!/bin/bash
multi_tag \
 --mapfile ./tag.mapfile\
 --mod shell\
 --thread 15 \
 --barcode_fasta ../Clindex/tag_barcode.fasta\
 --fq_pattern L25C15\
 --split_fastq \
 --split_matrix
