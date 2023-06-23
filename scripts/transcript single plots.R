###hnrnpl - transcript level analysis
tx2gene = "~/Desktop/rbp_bed/gencode.v40.tx2gene.csv"
metadata_chx <- fread("~/Documents/GitHub/splicing_chx/data/metadata_chx.csv")

files = unique(file.path("~/Documents/GitHub/splicing_chx/data/salmon_chx/",metadata_chx$sample,"quant.sf")) 
names(files) = unique(metadata_chx$sample)

txi.tx <- tximport(files, 
                   type="salmon", 
                   tx2gene=tx2gene,
                   ignoreTxVersion = TRUE,
                   ignoreAfterBar = TRUE,
                   txOut = TRUE)
TPM_transcripts = as.data.frame(txi.tx$abundance) %>% 
  tibble::rownames_to_column(.,var="transcript_id") %>%
  mutate(enstxp = gsub("\\..*", "", transcript_id)) %>%
  left_join(annotables::grch38_tx2gene) %>%
  left_join(annotables::grch38[,c(1,3)])

hnrnpl <- TPM_transcripts %>%
  filter(symbol == "HNRNPL") %>%
  pivot_longer(c(2:17)) %>%
  mutate(group = gsub("_[^_]+$", "", name))

for_single_gene_plot <- hnrnpl %>%
  filter(#enstxp == "ENST00000221419" | 
    enstxp == "ENST00000388749") 

stattest <- for_single_gene_plot %>%
  group_by(enstxp) %>%
  t_test(value ~ group) %>%
  add_significance() %>% 
  add_xy_position(x="group")
stattest

for_single_gene_plot %>%
  ggplot(aes(x = group, y = value, color = enstxp)) +
  geom_boxplot(alpha = 0.1) +
  geom_point(position = position_dodge2(width = 1)) +
  #stat_pvalue_manual(stattest, hide.ns = T) + #[c(1,2),]) +
  #scale_y_continuous(breaks = c(0,100,1000,5000,10000,20000,30000)) +
  scale_y_log10() +
  #facet_wrap(facets = vars(symbol)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90))
  
    #(transcript == "ENST00000634753" | transcript == "ENST00000634983")

