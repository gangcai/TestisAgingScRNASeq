### basic analysis pipeline for Seurat V5 ###

#record program start time
start_time <- Sys.time()

#loading seurat library
library(Seurat)
options(future.globals.maxSize = 20000 * 1024^2) # 20G memory, Notice: change the memory usage accordingly
options(Seurat.object.assay.version = "v5") # use Seurat V5 format

### parameter settings ####
project_id <- "Testis"
min.features <- 200
min.cells <- 3
pc.num_integration <- 60 # number of PCs for integration
pc.num <- 20         # number of PCs used to clustering
resolution <- 0.3    # higher values will generate more clusters, and vice versa, recommended values for experiments: 0.05, 0.1, 0.15, 0.2, 0.25, 0.3
nCount_RNA.max <- 50000    # Notice: change this value according to your data QC results
nCount_RNA.min <- 1000
max.features <- 7000       # Notice: change this value according to your data QC results
per.mt <- 10  # maximal mitochondrial percentage


########## loading the 10X matrix files ##################
samples=c("TO1","TO2","TO3","TY1","TY2","TY3")
raw.data.merged=""
i=0
dir_base <- "/home/db/private/XieLab/aging/mouse_testis_scRNASeq//singleron_6samples//clean/matrix/"
for(sample in samples){
	dir <- paste0(dir_base,sample,"_matrix/")
        obj.counts <- Read10X(data.dir = dir)
        colnames=colnames(obj.counts)
        colnames(obj.counts)=paste0(sample,"_",colnames)
        i=i+1
        if(i == 1){
                raw.data.merged=obj.counts
        }else{
                raw.data.merged=cbind(raw.data.merged,obj.counts)
        }
}

all.cells=colnames(raw.data.merged)
all.samples=sapply(all.cells,function(x){strsplit(as.character(x),"_")[[1]][1]})
all.samples=as.character(all.samples)
groups <- sapply(all.samples,function(x){
  if(x == "TO1" | x == "TO2" | x == "TO3"){
    return("Old")
  }
  if(x=="TY1" | x == "TY2" | x == "TY3"){
    return("Young")
  }
})


metadata = data.frame("groups"=groups,"samples"=all.samples,"cells"=all.cells,row.names=all.cells)
write.table(metadata,file="primary.meta.data.tsv",sep="\t",row.names=F,col.names=T)

################## Seurat basic analysis pipeline ##################
obj <- CreateSeuratObject(raw.data.merged, meta.data = metadata, project = project_id, min.cells = min.cells, min.features = min.features)

obj <- PercentageFeatureSet(obj, pattern = "^mt-", col.name = "percent.mt")
p0 <- VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0, group.by = "samples")
pdf(paste0(project_id, "_sctransform_QC_before.pdf"))
print(p0)
dev.off()

#data filtering
obj <- subset(obj, subset = nFeature_RNA > min.features & nFeature_RNA < max.features & percent.mt < per.mt & nCount_RNA < nCount_RNA.max &  nCount_RNA > nCount_RNA.min )


p0 <- VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0, group.by = "samples")
pdf(paste0(project_id, "_sctransform_QC_after.pdf"))
print(p0)
dev.off()

### split by samples
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$samples)

### normalization and scale, dimension reduction
obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj)
obj <- RunPCA(obj,npcs=pc.num_integration)

### intergration , slow step
obj <- IntegrateLayers(
  object = obj, method = HarmonyIntegration,npcs=pc.num_integration,
  orig.reduction = "pca", new.reduction = "harmony",
  verbose = FALSE
)

### clustering
obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:pc.num)
obj <- FindClusters(obj, resolution = resolution, cluster.name = "clusters")
obj <- RunUMAP(obj, reduction = "harmony", dims = 1:pc.num, reduction.name = "umap")

saveRDS(obj,paste0(project_id,"_SeuratV5.rds"))

p1 <- DimPlot(
  obj,
  reduction = "umap",
  group.by = c("samples","clusters"),
  combine = T,
  ncol=2
)

pdf(paste0(project_id,"_umap.samplesMerged.pdf"),width=9,height=4)
print(p1)
dev.off()


p2 <- DimPlot(
  obj,
  reduction = "umap",
  group.by=c("clusters"),
  split.by="samples",
  combine = FALSE,
  ncol=2
)

pdf(paste0(project_id,"_umap_samplesSplitted.pdf"),width=8,height=6)
print(p2)
dev.off()


p3 <- DimPlot(
  obj,
  reduction = "umap",
  group.by = c("clusters"),
  label=T
)

pdf(paste0(project_id,"_umap_showClusters.pdf"),width=5,height=4)
print(p3)
dev.off()


############# Find marker genes for each cluster ###############
obj <- JoinLayers(obj)
markers <- FindAllMarkers(obj,
                             only.pos = TRUE, 
                             assay="RNA",
                             slot="data",
                             min.pct = 0.25, logfc.threshold = 0.25)

write.table(markers,file=paste0(project_id,"_merged_markers.tsv"),sep="\t",quote=F,row.names=F)

############ output cell ids for each cluster ###################
clusters=obj$seurat_clusters
clusters=data.frame(clusters)
write.table(clusters,file="cell_cluster_ids.tsv",sep="\t",row.names=T,col.names=F,quote=F)

############ output meta data ###########
metadata <-  obj@meta.data
write.table(metadata,file="meta.data.tsv",sep="\t",row.names=F,col.names=T)

############# Print sesssion info ############
writeLines(capture.output(sessionInfo()), "sessionInfo.txt")

#check elapsed time
end_time <- Sys.time()
print(capture.output(end_time-start_time))

##end of the pipeline
