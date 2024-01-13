library(Seurat)
obj <- readRDS("Testis_SeuratV5.rds")
anno_df <- read.table("annotation/cluster_annotation.tsv",header=T,sep="\t")
meta_data <- obj@meta.data
cellType <- sapply(meta_data$seurat_clusters,function(x){
        x <- as.character(x)
        anno_df[anno_df$Cluster == x,]$CellType
})
meta_data$CellType <- unlist(cellType)
new_idents <- meta_data$CellType
names(new_idents) <- rownames(meta_data)
new_idents <- factor(new_idents)
Idents(obj) <- new_idents
obj@meta.data <- meta_data
write.table(meta_data,file="meta.data.annotated.tsv",sep="\t",row.names=F,col.names=T)
saveRDS(obj,"Testis_SeuratV5_annotated.rds")
print("completed")
