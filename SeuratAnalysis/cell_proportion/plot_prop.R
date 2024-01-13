library(dplyr)
library(ggplot2)
library(viridis)
library("cowplot")

data <- read.table("Testis_Cell_Proportion_Changes.tsv",header=T,sep="\t")

data$CellType <- as.character(data$CellType)
cell_types <- unique(data$CellType)

groups <- unique(data$group_id)

results <- list("group"=c(),
                "cluster"=c(),
                "mean"=c(),
                "sd"=c(),
                "celltype"=c())

fcs <- list("celltype"=c(),
            "logfc"=c())
for(cell_type in cell_types){
  data_s <- data[data$CellType == cell_type, ]
  
  youngs <- data_s[data_s$group_id == "TY",]$cell_percentage
  olds <- data_s[data_s$group_id == "TO",]$cell_percentage
  young_mean <- mean(youngs)
  old_mean <- mean(olds)
  young_sd <- sd(youngs)
  old_sd <- sd(olds)
  results[["group"]] <- c(results[["group"]] , c("young","old"))
  results[["cluster"]] <- c(results[["cluster"]],c(cell_type,cell_type))
  results[["mean"]] <- c(results[["mean"]],c(young_mean,old_mean))
  results[["sd"]] <- c(results[["sd"]],c(young_sd,old_sd))
  results[["celltype"]] <- c(results[["celltype"]],c(cell_type,cell_type))
  logfc <- log2(old_mean/young_mean)
  fcs[["celltype"]] <- c(fcs[["celltype"]],cell_type)
  fcs[["logfc"]] <- c(fcs[["logfc"]],logfc)
  
}



fcs_df <- bind_rows(fcs)
fcs_df <- arrange(fcs_df, logfc)
fcs_df$celltype <- factor(fcs_df$celltype, levels=fcs_df$celltype)
p1 <- ggplot(fcs_df, aes(celltype , logfc, fill=celltype)) +
   geom_bar(stat="identity") +
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1),
        legend.position="none")+
  ylab("log2(percentage foldchange)") +
  xlab("") +
  scale_fill_viridis(discrete = TRUE) 


pdf("CellProportion_FC_barplot.pdf",width=9,height=5)
print(p1)
dev.off()


results_df <- bind_rows(results)
results_df$celltype <- factor(results_df$celltype, levels=fcs_df$celltype)
p2 <- ggplot(results_df, aes(x=celltype, y=mean, fill=group)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,
                position=position_dodge(.9)) +
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1))+
  xlab("")+
  ylab("Cell Proportion")
  
pdf("CellProportion_Barplot.pdf",width=9,height=5)
print(p2)
dev.off()

p <- plot_grid(p2,p1,ncol=1)
pdf("CellProportion_Barplot_FC_combined.pdf")
print(p)
dev.off()
