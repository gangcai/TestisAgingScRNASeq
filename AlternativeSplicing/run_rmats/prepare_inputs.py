import re
import glob
celltypes=[]
for filename in glob.glob("../bam/TO1/*.bam"):
	celltype = filename.split("/")[-1]
	celltype = celltype.split(".")[1]
	celltypes.append(celltype)

print(celltypes)
for celltype in celltypes:
	filenames = []
	for sample in ["TO1","TO2","TO3"]:
		filename="../bam/" + sample + "/" + sample + "." + celltype + ".bam"
		filenames.append(filename)	
	outfile = open(celltype + "_old.txt","w")
	outfile.write(",".join(filenames)+"\n")
	outfile.close()
	filenames = []
	for sample in ["TY1","TY2","TY3"]:
		filename="../bam/" + sample + "/" + sample + "." + celltype + ".bam"
		filenames.append(filename)	
	outfile = open(celltype + "_young.txt","w")
	outfile.write(",".join(filenames)+"\n")
	outfile.close()
	
