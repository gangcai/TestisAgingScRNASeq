import re
i = 0
sample2info={}
for line in open("../../basic_analysis/meta.data.annotated.tsv"):
	i += 1
	if i==1:
		continue
	line=re.sub("\"","",line)
	items = line.rstrip().split("\t")
	cells_with_sample = items[3]
	(sample,cell_name)=cells_with_sample.split("_")
	celltype = items[-1]
	if sample in sample2info:
		sample2info[sample].append(cell_name+"\t"+celltype)
	else:
		sample2info[sample]=[cell_name+"\t"+celltype]

for sample in sample2info:
	outfile=open(sample+"_cell_cluster_info.tsv","w")
	outfile.write("Index\tCell_type\n")
	infos = sample2info[sample]
	for info in infos:
		outfile.write(info+"\n")
outfile.close()

