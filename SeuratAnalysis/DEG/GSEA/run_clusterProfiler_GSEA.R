library(AnnotationHub)
library(clusterProfiler)
library(ReactomePA)
library(MeSHDbi)
library(meshes)
deg_df <- read.table("../../TestisAge_DEGs_WilCox.tsv",header=T,sep="\t")
deg_df$cellType <- deg_df$cluster.info

for(cid in unique(deg_df$cellType)){
  deg_df_s <- deg_df[deg_df$cellType == cid,]
  
  deg_df_ordered <- deg_df_s[order(deg_df_s$avg_log2FC,decreasing=T),]
  
  genes_ranked <- deg_df_ordered$gene
  genes_ranked_ids <- bitr(genes_ranked, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
  
  fc_ranked <- deg_df_ordered$avg_log2FC
  
  names(fc_ranked) <- genes_ranked
  
  #run gseGO
  go_results <- gseGO(geneList     = fc_ranked,
                      OrgDb        = "org.Mm.eg.db",
                      ont          = "ALL",
                      minGSSize    = 30,
                      maxGSSize    = 500,
                      pvalueCutoff = 0.1,
                      keyType="SYMBOL",
                      verbose      = FALSE)
  
  write.table(data.frame(go_results),file=paste0(cid,"_GSEA_GO.tsv"),sep="\t",row.names = F,col.names = T)
  
  #run GSEA KEGG
  fc_ranked_2 <- fc_ranked[genes_ranked_ids$SYMBOL]
  names(fc_ranked_2) <- genes_ranked_ids$ENTREZID
  
  kegg_results <- gseKEGG(fc_ranked_2, 
                          organism = "mmu",
                          minGSSize = 10,
                          pvalueCutoff=0.1,
                          keyType="ncbi-geneid")
  
  kegg_results_readable <- setReadable(kegg_results, 'org.Mm.eg.db', keyType = "ENTREZID")
  
  write.table(kegg_results_readable@result,file=paste0(cid,"_GSEA_KEGG.tsv"),sep="\t",row.names = F,col.names = T)
  
  #run reactome
  reactome_results <- gsePathway(fc_ranked_2, 
                                 pvalueCutoff = 0.1,
                                 organism="mouse",
                                 pAdjustMethod = "BH", 
                                 verbose = FALSE)
  reactome_results_readable <- setReadable(reactome_results, 'org.Mm.eg.db', keyType = "ENTREZID")
  
  write.table(reactome_results_readable@result,file=paste0(cid,"_GSEA_Reactome.tsv"),sep="\t",row.names = F,col.names = T)
  
}
