library(dplyr)
library(ggplot2)
p_adj_cutoff <- 0.05
logfc_cutoff <- log2(2)
percentage_cutoff <- 0.1

project_id <- "TestisAge"
data <- read.table(paste0("../",project_id,"_DEGs_WilCox_BH.tsv"),header=T,sep="\t")
data.f <- data %>% filter(p_val_adj_BH < p_adj_cutoff & abs(avg_log2FC) > logfc_cutoff  & (pct.1 > percentage_cutoff | pct.2 > percentage_cutoff))
write.table(data.f,file=paste0(project_id,"_DEGs_WilCox_Sig.tsv"),col.names=T,row.names=F,sep="\t",quote=F)

#data.f$cluster.info <- factor(as.numeric(data.f$cluster.info))
stats1 <- data.f %>% 
    group_by(cluster.info) %>%
    summarise(DEG_Num=n()) %>% arrange(desc(DEG_Num))

colnames(stats1) <- c("cluster","DEG_Num")
write.table(stats1,file=paste0(project_id,"_DEGs_WilCox_Sig_Summary.tsv"),
	    col.names=T,row.names=F,sep="\t",quote=F)

celltype_o <- stats1$cluster

stats_up <- data.f %>% filter(avg_log2FC > 0) %>% group_by(cluster.info) %>% summarise(DEGNum=n())

stats_down <- data.f %>% filter(avg_log2FC < 0) %>% group_by(cluster.info) %>% summarise(DEGNum=n())

stats_up$changeType <- "UP"
stats_down$changeType <- "DOWN"

stats <- rbind(stats_up,stats_down)

stats$cluster.info <- factor(stats$cluster.info,levels=celltype_o)
p <- ggplot(stats, aes(x=cluster.info, y=DEGNum, fill=changeType)) + geom_bar(stat="identity") +
  theme(axis.text.x = element_text(angle = 60, vjust = 0.8, hjust=0.8)) +
  xlab("")+ylab("Number of age-DEGs")+
  scale_fill_manual(values=c("darkseagreen","coral"))

write.table(stats,file=paste0(project_id,"_DEGs_WilCox_Sig_Summary_UpDown.tsv"),
            col.names=T,row.names=F,sep="\t",quote=F)

pdf("TestisAge_DEG_Summary.pdf",width=6,height=4)
print(p)
dev.off()


