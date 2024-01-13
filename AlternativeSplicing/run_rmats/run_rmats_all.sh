#!/bin/bash
gtf="/home/db/public/SoftwareIndex/celescope1.16.1/mouse/gencode.vM25.primary_assembly.annotation.gtf"
#celltype="Early_diffSPG"
for filename in *_old.txt
do
	celltype=${filename//_old.txt/}
	echo ${celltype}
	b1=${celltype}_old.txt
	b2=${celltype}_young.txt
	outdir=${celltype}
	mkdir -p $outdir
	outdir_tmp=${outdir}.tmp
	mkdir -p ${outdir_tmp}
	echo "begin"
	date
	rmats --b1 $b1 --b2 $b2 --gtf $gtf --od $outdir --tmp $outdir_tmp -t single --readLength 150 --nthread 20 --tstat 10
	echo "completed"
	date
done
