library(Seurat)
library(dplyr)
library(ggplot2)
obj <- readRDS("../Testis_SeuratV5_annotated.rds")

meta.data <- obj@meta.data
meta.data$cell <- rownames(meta.data)


stats1 <- meta.data %>% group_by(samples) %>% summarise(cellNum=n())

stats2 <- meta.data %>% group_by(CellType) %>% summarise(cellNum=n()) %>% arrange(desc(cellNum))

stats3 <- meta.data %>% group_by(samples,CellType) %>% summarise(cellNum=n())

write.table(stats1,file="TestisAging_CellNumberStats_Samples.tsv",
            row.names = F,
            sep="\t")

write.table(stats2,file="TestisAging_CellNumberStats_CellType.tsv",
            row.names = F,
            sep="\t")

write.table(stats3,file="TestisAging_CellNumberStats_CellTypeSamples.tsv",
            row.names = F,
            sep="\t")



celltypes_o <- stats2$CellType

stats3$CellType <- factor(stats3$CellType,levels = celltypes_o)

p <- ggplot(stats3, aes(CellType,cellNum,fill=CellType)) + geom_bar(stat = "identity") + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1), legend.position = "none")+
  labs(title = NULL, x = NULL, y = "Number of cells") +
  facet_wrap(~samples)

pdf("Barplot_Samples_CellType.pdf",width = 8,height=4)
print(p)
dev.off()


p <- ggplot(stats1,aes(samples,cellNum,fill=samples)) + geom_bar(stat = "identity") +
  theme(legend.position = "none")+
  labs(title = NULL, x = NULL, y = "Number of cells")

pdf("Barplot_Samples.pdf",width = 4,height=4)
print(p)
dev.off()