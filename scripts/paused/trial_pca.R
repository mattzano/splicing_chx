### pca analysis
pca <- normed_counts_chx[,c(1:9)] %>%
  arrange(CTRL_chx_1) %>%
  distinct(Geneid, .keep_all = TRUE) %>%
  #select(c(1:)) %>%
  column_to_rownames('Geneid')
#column_to_rownames('ensgene') #%>% t()
#pca_res <- prcomp(t(pca),scale = F)

meta_df <- metadata_dir %>% 
  filter(treatment == "chx") %>%
  #mutate(condition = factor(condition, levels = c(0,1))) %>%
  #mutate(treatment = factor(treatment, levels = c(0,1))) %>%
  mutate(group = gsub("_[^_]+$", "", sample)) %>%
  column_to_rownames('sample') 

p <- PCAtools::pca(pca, metadata = meta_df, removeVar = 0.1)

screeplot(p, axisLabSize = 18, titleLabSize = 22)

#PCAtools::eigencorplot(p,components = getComponents(p),
#                       metavars = 'condition')

biplot(p,
       x = "PC1",
       y = "PC2",
       colby = 'group')

pc1_gene_loadings = p$loadings %>% 
  dplyr::select(PC1) %>% 
  rownames_to_column('Geneid') %>% 
  #left_join(annotables::grch38 %>% dplyr::select(ensgene)) %>% 
  #https://statisticsglobe.com/r-dplyr-join-inner-left-right-full-semi-anti
  arrange(-PC1) 

for_ggplot <- normed_counts_chx[,c(1:9,11)] %>% 
  left_join(pc1_gene_loadings %>% head(16)) %>%
  filter(!is.na(PC1)) %>%
  #dplyr::select(-ensgene) %>% 
  melt(id.vars = c("Geneid", "symbol","PC1")) %>%
  mutate(Geneid = fct_reorder(Geneid,-PC1)) %>% 
  separate(variable, into = c("condition","treatment"))

for_ggplot %>%
  ggplot(aes(x = treatment, y = value, fill = condition)) + 
  geom_col(position = "dodge2") + 
  facet_wrap(~symbol)


##this tool is not working - would be good cause it includes lot of data and GO
#annotation <- normed_counts[,c(1,19)] %>%
#  distinct(Geneid, .keep_all = TRUE) %>%
#  column_to_rownames('Geneid')
#se <- DESeqDataSetFromMatrix(countData = pca,
#                              colData = meta_df,
#                              design = ~ condition + treatment)
#se <- DESeqTransform(se)
#go <- pcaExplorer::pca2go(
#    se,
#  annotation = annotation,
#  inputType = "geneSymbol",
#  organism = "Hs",
#  ensToGeneSymbol = TRUE
#)
#pcaExplorer(countmatrix = pca, coldata = meta_df, annotation = annotation, pca2go = go)
