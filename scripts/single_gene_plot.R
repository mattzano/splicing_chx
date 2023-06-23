for_single_gene_plot <- normed_counts_UPF1_long %>%
  dplyr::filter(symbol == "TARDBP") # | symbol == "UNC13A" | symbol == "STMN2" | symbol == "AARS")

library(rstatix)
stattest <- for_single_gene_plot %>%
  group_by(symbol) %>%
  t_test(value ~ group) %>%
  add_significance() %>% 
  add_xy_position(x="group")
stattest

for_single_gene_plot %>%
  ggplot(aes(x = group, y = value, color = symbol)) +
  geom_boxplot(alpha = 0.1) +
  geom_point(position = position_dodge2(width = 1)) +
  #stat_pvalue_manual(stattest, hide.ns = T) + #[c(1,2),]) +
  #scale_y_continuous(breaks = c(0,100,1000,5000,10000,20000,30000)) +
  scale_y_log10() +
  facet_wrap(facets = vars(symbol)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90))



#compare chx and upf1
AppendMe <- function(dfNames) {
  do.call(rbind, lapply(dfNames, function(x) {
    cbind(get(x), source = x)
  }))
}

normed_long_both <- AppendMe(c("normed_counts_chx_long", "normed_counts_UPF1_long")) 

normed_long_both %>%
  dplyr::filter(symbol == "TARDBP" | symbol == "STMN2" | symbol == "AARS" | symbol == "UNC13A" | symbol == "UPF2") %>%
  ggplot(aes(x = group, y = value, color = source)) +
  geom_boxplot(alpha = 0.1, show.legend = F) +
  geom_point(position = position_dodge2(width = 1), show.legend = F) +
  #stat_pvalue_manual(stattest, hide.ns = T) + #[c(1,2),]) +
  #scale_y_log10() +
  facet_wrap(facets = vars(symbol), scales = "free_y") +
  scale_color_brewer(palette = "Set1") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90))

###see if more of a general trend
results_table_ctrl_both <- full_join(results_table_cont, results_table_UPF1_ctrl, by = "symbol", na_matches = "never")
results_table_nmd_both <- full_join(results_table_cyhx, results_table_UPF1_UPF1KD, by = "symbol")
#results_table_cont_norm <- results_table_cont %>%
#  mutate(normalised_log2_chx = (log2FoldChange - mean(log2FoldChange)) / sd(log2FoldChange))
#results_table_UPF1_ctrl_norm <- results_table_UPF1_ctrl %>%
#  mutate(normalised_log2_upf1 = (log2FoldChange - mean(log2FoldChange)) / sd(log2FoldChange))
#results_table_ctrl_both <- left_join(results_table_cont_norm, results_table_UPF1_ctrl_norm, by = "symbol")

#results_table_ctrl_both %>%
results_table_ctrl_both %>%
  dplyr::filter(`padj.x` < 0.05 & `padj.y` < 0.05) %>%
  #mutate(normalised_log2_chx = (abs(log2FoldChange.x) - min(abs(log2FoldChange.x)) / (max(abs(log2FoldChange.x) - min(abs(log2FoldChange.x)))))) %>%
  #mutate(normalised_log2_upf1 = (abs(log2FoldChange.y) - min(abs(log2FoldChange.y)) / (max(abs(log2FoldChange.y) - min(abs(log2FoldChange.y)))))) #%>%
  ggplot(aes(x = `log2FoldChange.x`, y = `log2FoldChange.y`)) +
  geom_point() +
  geom_smooth(method = 'lm', se = F) +
  stat_cor(method="pearson") +
  scale_color_brewer(palette = "Set1") +
  xlab('Log2 Fold Change in CHX experiment') +
  ylab('Log2 Fold Change in UPF1KD experiment') +
  #coord_fixed() +
  theme_classic()

ggsave("results/upgrade/correlation_expression_ctrl.png")
ggsave("results/upgrade/correlation_expression_nmd.png")

ds_both <- full_join(`DS_ControlControl-ControlTDP43KD`, `DS_CtrlCtrl-TDP43Ctrl`, by = "paste_into_igv_junction") 
ds_both <- full_join(`DS_CycloheximideControl-CycloheximideTDP43KD`, `DS_CtrlUPF1-TDP43UPF1`, by = "paste_into_igv_junction") 


list_gene <- c("UNC13A", "STMN2") #, "AARS1", 
               #"CELF5", 
               #"ELAVL3", "INSR", #"UNC13B", 
               #"ATG4B", "KCNQ2", "SYT7", "WASL", "RAP1GAP", "AGRN","PFKP", 
               #"HDGFL2","ARHGAP32", "CEP290", "PRELID3A", "KALRN",
               #"CAMK2B", "ADCY1", "TRAPPC12", "EPB41L4A", "IQCE", "SEMA4D"
               #"CYFIP2", "SYNE1")

splicing_corr <- ds_both %>%
  dplyr::filter((#`probability_changing.x` > 0.9 & 
                 `base_mean_psi.x` < 0.05 & `mean_dpsi_per_lsv_junction.x` > 0.1) &
                  (#`probability_changing.y` > 0.9 & 
                   `base_mean_psi.y` < 0.05 & `mean_dpsi_per_lsv_junction.y` > 0.1)) %>%
  #mmutate(`Novel Junction` = de_novo_junctions == 1) %>%
  #mutate(log10_test_stat = -log10(1 - probability_changing)) %>%
  #mutate(log10_test_stat = ifelse(is.infinite(log10_test_stat), 4.5, log10_test_stat)) %>%
  #mutate(graph_alpha = ifelse(probability_changing > 0.9, 1, 0.2)) %>%
  mutate(label_junction = case_when(`gene_name.x` %in% list_gene &
                                      `de_novo_junctions.x` == 1 ~ `gene_name.x`,
                                    T ~ "")) 

#788(314) 533(320) 183(134) 1138(500) 
splicing_corr %>%
  ggplot(aes(x = `mean_dpsi_per_lsv_junction.x`, y = `mean_dpsi_per_lsv_junction.y`)) +
  geom_point() +
  geom_smooth(method = 'lm', se = F) +
  geom_text_repel(aes(label = label_junction),
                  point.padding = 0.3,
                  nudge_y = -0.2,
                  min.segment.length = 0,
                  box.padding  = 3,
                  max.overlaps = Inf,
                  size=4, show.legend = F) +
  stat_cor(method="pearson") +
  scale_color_brewer(palette = "Set1") +
  coord_fixed() +
  xlab(expression(paste(Delta, Psi,' after TDP-43 KD in SH-SY5Y'))) +
  ylab(expression(paste(Delta, Psi,' after TDP-43 KD in i3Neurons'))) +
  theme_classic() +
  theme(text = element_text(size = 22))

ggsave(filename ="~/Desktop/correlation_splicing_nmd.png", width = 7)
