library(Seurat)
library(ggplot2)
library(dplyr)


deg_df <- read.table("../summary_BH/TestisAge_DEGs_WilCox_Sig.tsv",header=T,sep="\t")


obj <- readRDS("../../../Testis_SeuratV5_annotated.rds")


for(celltype in c("Telocyte","EC","myoid","Leydig","macrophage","Tcell")){
  top_up_degs <- deg_df %>% filter(cluster.info == celltype & avg_log2FC >0) %>% top_n(5,-p_val_adj_BH) %>% top_n(5,avg_log2FC)
  top_down_degs <- deg_df %>% filter(cluster.info == celltype & avg_log2FC <0) %>% top_n(5,-p_val_adj_BH) %>% top_n(5,-avg_log2FC)
  
  top_degs <- rbind(top_up_degs,top_down_degs)
  
  
  p <- VlnPlot(obj, 
               features = top_degs$genes,
               split.by = "groups",
               idents = celltype,
               ncol=5)
  pdf(paste0(celltype,"_topDEGs.pdf"),width = 10,height=5)
  print(p)
  dev.off()
  
  
  p <- VlnPlot(obj, 
               features = top_degs$genes,
               split.by = "groups",
               ncol=2)
  pdf(paste0(celltype,"_topDEGs_showAllCellTypes.pdf"),width = 10,height=15)
  print(p)
  dev.off()
}
