for prefix in ["TY1"]:
	sample="TestisMix"
	fpath="../transcriptome_seq/P22033104/rawdata/"
	outfile=open(sample+".mapfile","w")
	outfile.write("\t".join([prefix,fpath,sample])+"\n")
	outfile.close()
