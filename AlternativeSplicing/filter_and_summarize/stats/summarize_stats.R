library(dplyr)

as_types <- c("SE","A3SS","A5SS","MXE","RI")

all_data <- data.frame()
for(as_type in as_types){
  filename <- paste0("../",as_type,"_significant_summary.tsv")
  data <- read.table(filename,sep="\t",header=T)
  
  data$AS_Type <- as_type
  
  if(as_type == "SE"){
    all_data <- data
  }else(
    all_data <- rbind(all_data,data)
  )
}

write.table(all_data,file="TestisAging_Sig_AS_summary.tsv",
            sep="\t",
            row.names = F,
            quote=F)

stats1 <- all_data %>% group_by(AS_Type) %>% summarise(counts=sum(SigEventsCounts))


### get the number of unique genes
as_genes_data <- read.table("../GeneEnrichmentAnalysis/All_Sig_AS_Genes.tsv",
                            header=T,
                            sep="\t")

unique_as_genes <- unique(as_genes_data$gene)

unique_df <- data.frame(gene=unique_as_genes)
write.table(unique_df,file="unique_AS_genes.tsv",row.names = F,sep="\t",quote=F)



