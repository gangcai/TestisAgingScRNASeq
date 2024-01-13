library(Seurat)
library(cowplot)
library(dplyr)
obj = readRDS("../Testis_SeuratV5_annotated.rds")

p <- DimPlot(obj, reduction = "umap",pt.size =0.1, group.by = "samples")
pdf("Testis_UMAP_merged_SampleColored.pdf",width=8,height=7)
print(p)
dev.off()

p <- DimPlot(obj, reduction = "umap",pt.size =0.1, split.by = "samples",ncol=2) + NoLegend()
pdf("Testis_UMAP_merged_SampleSplitted.pdf",width=4,height=4)
print(p)
dev.off()

p <- DimPlot(obj, label = TRUE, reduction="umap") + NoLegend()
pdf("Testis_UMAP_merged.pdf",width=7,height=7)
print(p)
dev.off()

p <- DimPlot(obj, label = F, reduction="umap")
pdf("Testis_UMAP_merged_legend.pdf",width=8,height=7)
print(p)
dev.off()
