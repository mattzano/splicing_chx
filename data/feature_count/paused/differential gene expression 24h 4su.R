results <- my.dds$results_table
filtered_res <- filter(results, padj < 0.1 )
ordered_res <- arrange(filtered_res, desc(abs(log2FoldChange)))
#top_30 <- head(ordered_res, 30) %>% 
#  select(gene_name)
#add_name_to_plot(my.dds$results_table, gene_names = top_30$gene_name)

#Now I am doing a GO analysis on this data

library(clusterProfiler)
library(org.Hs.eg.db)

gene <- ordered_res$gene_name[abs(ordered_res$log2FoldChange) > 2]

ggo <- groupGO(gene = gene,
               OrgDb = org.Hs.eg.db,
               keyType = "SYMBOL",
               ont = "CC",
               level = 3)
head (ggo)

ego <- enrichGO(gene = gene,
                keyType = "SYMBOL",
                OrgDb = org.Hs.eg.db,
                ont = "CC",
                pAdjustMethod = "BH",
                pvalueCutoff = 0.01,
                qvalueCutoff = 0.05)
ego_df <- ego@result
ego_df_sort <- arrange(ego_df, p.adjust)

#gene.df <- bitr(gene, fromType = "SYMBOL",
#                toType = c("ENSEMBL","ENTREZID"),
#                OrgDb = org.Hs.eg.db, drop = T)

#ego2 <- enrichGO(gene = gene.df$ENSEMBL,
#                 OrgDb = org.Hs.eg.db,
#                 keyType = 'ENSEMBL',
#                 ont = "CC",
#                 pAdjustMethod = "BH",
#                 pvalueCutoff = 0.01,
#                 qvalueCutoff = 0.05)
#head(ego2, 3)

#ego3 <- gseGO(geneList = gene,
#              OrgDb = org.Hs.eg.db,
#              ont = "CC",
#              minGSSize = 100,
#              maxGSSize = 500,
#              pvalueCutoff = 0.05,
#              verbose = FALSE)

goplot(ego)
