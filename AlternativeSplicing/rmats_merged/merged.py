import glob,re
as_types=["SE","RI","MXE","A5SS","A3SS"]
outs={}
records={}	
for as_type in as_types:
	outs[as_type]=open(as_type+"_merged.tsv","w")
	records[as_type]=0

for filename in glob.glob("../run_rmats/*_old.txt"):
	sample=filename.split("/")[-1]
	sample=re.sub("_old.txt","",sample)
	for rf in glob.glob("../run_rmats/"+sample+"/*MATS.JCEC.txt"):
		as_type=rf.split("/")[-1]
		as_type=re.sub(".MATS.JCEC.txt","",as_type)
		records[as_type]+=1
		i=0
		for line in open(rf):
			i+=1	
			if i==1 and records[as_type]==1:
				outs[as_type].write("cellType\t"+line)
				continue
			if i==1:
				continue
			outs[as_type].write(sample+"\t"+line)
for as_type in as_types:
	outs[as_type].close()
