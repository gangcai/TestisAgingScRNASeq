library(dplyr)
library(ggplot2)

for(as_type in c("SE","A3SS","A5SS","MXE","RI")){
  data <- read.table(paste0("../",as_type,"_merged.tsv"),header=T,sep="\t")
  
  data_sig <- data %>% filter((FDR < 0.05) & abs(IncLevelDifference) > 0.2) %>%
    arrange(FDR)
  
  write.table(data_sig,file=paste0(as_type,"_significant.tsv"),row.names = F,col.names=T,sep="\t")
  
  data_summary <- data_sig %>% group_by(cellType) %>% summarise(SigEventsCounts=n()) %>%
    arrange(desc(SigEventsCounts))
  
  write.table(data_summary,file=paste0(as_type,"_significant_summary.tsv"),row.names = F,col.names=T,sep="\t")
  
  data_summary$cellType <- factor(data_summary$cellType, levels=data_summary$cellType)
  p <- ggplot(data_summary, aes(cellType, SigEventsCounts)) + geom_bar(stat="identity",
                                                                       position=position_dodge(),
                                                                       fill="lightpink") +
    theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) +
    xlab("") +
    ylab("# of significant AS events")
  
  pdf(paste0(as_type, "_significant_summary_barplot.pdf"),width=5,height=4)
  print(p)
  dev.off()
  
  
  
}
