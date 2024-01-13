library(dplyr)
library(ggplot2)
meta_data <- read.table("../meta.data.annotated.tsv",header=T,sep="\t")

CellTypes <- unique(meta_data$CellType)
CellTypes <- as.character(CellTypes)

samples <- unique(meta_data$samples)
groups <- unique(meta_data$group)


#get number of cells for each sample
sample_cellnum_list <- table(meta_data$samples)

######### get number of cells for each sample in each CellType ############
cell_num_list <- list()

for(cid in CellTypes){
  data_f <- sapply(samples, function(x){
    meta_data_f <- meta_data %>%
      filter(CellType == cid, samples == x)
    cells_f <- meta_data_f$cell
    cell_num <- length(cells_f)
    
    sample_total_cell <- sample_cellnum_list[[x]]
    
    cell_per <- cell_num/sample_total_cell
    
    
    group_id <- gsub("[0-9]","",x)
    group_id <- gsub("-M","",group_id)
    result <- c(cid, group_id, cell_num,cell_per)
    names(result) <- c("CellType_id","group_id","cell_num","cell_percentage")
    return(result)
  })
  data_f <- data.frame(t(data_f))
  cell_num_list[[cid]] <- data_f
}

######### t-test ###########
for(cid in CellTypes){
  
   df <- cell_num_list[[cid]]
   num_1 <- as.numeric(df[df$group_id == "TO",]$cell_percentage)
   num_2 <- as.numeric(df[df$group_id == "TY",]$cell_percentage)
   fc <- log2(mean(num_1)/mean(num_2))
   p <- t.test(num_1,num_2)$p.value
   cell_num_list[[cid]]$log2_FC <- fc
   cell_num_list[[cid]]$p_value <- p
}



cell_num_tbl <- bind_rows(cell_num_list)
write.table(cell_num_tbl, file = "Testis_Cell_Proportion_Changes.tsv",sep="\t",row.names = F, col.names=T, quote=F)


######## Plot the results #############
# reload data 
#cell_num_tbl <- read.table("Testis_Cell_Proportion_Changes.tsv",header=T,sep="\t")
CellTypes <- unique(cell_num_tbl$CellType_id)
CellTypes <- as.character(CellTypes)
CellTypes <- CellTypes[order(CellTypes)]
cell_num_tbl$cell_percentage <- 100*as.numeric(cell_num_tbl$cell_percentage)
cell_num_tbl$CellType_id <- as.character(cell_num_tbl$CellType_id)
cell_num_tbl$CellType_id <- factor(cell_num_tbl$CellType_id, levels=as.character(CellTypes))

p <- ggplot(cell_num_tbl, aes(x = CellType_id, y = cell_percentage, fill=group_id)) +
     geom_boxplot() +
     xlab("Cell Types") +
     ylab("Cell Percentage")
pdf("Testis_Cell_Proprotion_Changes_Boxplot.pdf", width=10,height=4)
print(p)
dev.off()
