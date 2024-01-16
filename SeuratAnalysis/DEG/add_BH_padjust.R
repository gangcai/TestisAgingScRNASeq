data <- read.table("TestisAge_DEGs_WilCox.tsv",header = T,sep="\t")
data$p_val_adj_BH <- p.adjust(data$p_val,method="BH")
write.table(data,file="TestisAge_DEGs_WilCox_BH.tsv",
            row.names = F,
            col.names = T,
            quote=F,
            sep="\t")