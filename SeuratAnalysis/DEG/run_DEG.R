library(Seurat)
library(ggplot2)
######## cell-level wilcoxon two group DEG comparison ########
### settings ####
min_cell <- 10 # minimal number of cells for each group
project_id <- "TestisAge" # change this accroding to your project name

#### loading the data ####
obj <- readRDS("../../Testis_SeuratV5_annotated.rds") # change to your own rds file
DefaultAssay(obj) <- "RNA"

#### add group information #####
Idents(obj) <- "groups"
clusters <- unique(obj@meta.data$CellType)
groups <- c("Old","Young") # change this according to your group information

#### cluster-level wilcoxon test ######
obj <- JoinLayers(obj)
results <- ""
cells <- names(obj$groups)
i <- 0
for(cid in clusters){
   print(cid)
   group1 <- groups[1]
   cell.g1 <- cells[obj$groups==group1 & obj$CellType == cid]
   group2 <- groups[2]
   cell.g2 <- cells[obj$groups == group2 & obj$CellType == cid]
   if(length(cell.g1) > min_cell & length(cell.g2) > min_cell){
	   markers <- FindMarkers(obj,ident.1=cell.g1, # ident.1 vs ident.2
				  ident.2=cell.g2,
				  slot="data",
				  verbose=F,
			          assay="RNA",
				  test.use="wilcox",
				  logfc.threshold=0,
				  only.pos=F,
				  pseudocount.use=0.1)
	   cluster.info <- rep(cid,nrow(markers))
	   cmp_sample <- paste0(group1,"_vs_",group2)
	   cmp_info <- rep(cmp_sample,nrow(markers))
	   hgr <- cbind(cluster.info,markers)
	   hgr <- cbind(cmp_info,hgr)
	   genes <- rownames(hgr)
	   hgr <- cbind(genes,hgr)
	   i <- i+1
	   if(i == 1){
	     results <- hgr
	   }else{
	     results <- rbind(results,hgr)
	   }
   }
}

##### generate the output file #######
write.table(results,
	    file=paste0(project_id,"_DEGs_WilCox.tsv"),
	    sep="\t",
	    quote=F,
	    row.names=F,
	    col.names=T)
