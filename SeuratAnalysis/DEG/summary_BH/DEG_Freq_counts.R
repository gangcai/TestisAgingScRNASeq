library(ggplot2)
library(dplyr)

data <- read.table("TestisAge_DEGs_WilCox_Sig.tsv",
                   header=T,
                   sep = "\t")

gene_freq <- table(data$genes)

freq_df <- data.frame(gene=names(gene_freq),
                      count=as.numeric(gene_freq))

freq_df <- freq_df %>% arrange(desc(count))

write.table(freq_df,
            file="age_DEGs_counts.tsv",
            sep="\t",
            row.names = F,
            quote = F)