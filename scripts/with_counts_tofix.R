big_data_ce <- big_data[,c(1,2,5,17)] %>%
  pivot_wider(names_from = .id, values_from = mean_psi_per_lsv_junction, values_fn = mean) %>% ###239608
  dplyr::filter(Control_Control < 0.05 & Cycloheximide_Control < 0.05) %>% ##95487
  dplyr::filter(Control_TDP43KD > 0.1 | Cycloheximide_TDP43KD > 0.1) %>% ##2089
  #dplyr::filter(Cycloheximide_TDP43KD < 0.05) %>%
  dplyr::mutate(code = paste0(gene_name, "_", paste_into_igv_junction))

big_data_ce %>%
  dplyr::group_by(gene_name, paste_into_igv_junction) %>%
  dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
  dplyr::filter(n > 1L) 


big_data_ce_plot <- big_data_ce %>%
  #drop_na() %>%
  #rownames_to_column("code") %>%
  pivot_longer(cols = c(Control_Control,
                        Cycloheximide_Control, 
                        Control_TDP43KD, 
                        Cycloheximide_TDP43KD))

big_data_ce_plot$name <- factor(big_data_ce_plot$name, 
                                levels = c("Control_Control",
                                                         "Cycloheximide_Control", 
                                                         "Control_TDP43KD", 
                                                         "Cycloheximide_TDP43KD"))


plot <- big_data_ce_plot %>%
  #dplyr::filter(value < 2 & value > -2) %>%
  #dplyr::filter(paste_into_igv_junction == "chr8:79611214-79616822") %>%
  ggplot(mapping = aes(x = name, y = value, group = paste_into_igv_junction)) +
  geom_line(aes(), show.legend = F) +
  #geom_bar(fill ="#3288BD", stat = "identity", width = 0.5) +
  #geom_errorbar(aes(x = .id, ymin = mean_psi_per_lsv_junction - (stdev_psi_per_lsv_junction/sqrt(2)),
  #                           ymax = mean_psi_per_lsv_junction + (stdev_psi_per_lsv_junction/sqrt(2))), width = 0.3) +
  #geom_text_repel(aes(x = .id, y = mean_psi_per_lsv_junction, label = label_junction), point.padding = 0.3,
  #                nudge_y = 0.2, min.segment.length = 0.5, box.padding  = 2, max.overlaps = Inf, size=4, show.legend = F) +
  #scale_color_manual(values = c("gray","#3288BD",
  #                              "#D53E4F")) +
  #scale_alpha_manual(values = c(0.1,1)) +
  #xlab("TDP-43 protein levels") +
  ylab("PSI") +
  #scale_x_discrete(labels = c("100%", "77%", "25%", "8%", "4%", "0%")) +
  scale_y_continuous(breaks=seq(0,1,0.2),
                     labels = c("0%","20%","40%","60%","80%","100%")) +
  theme_classic()
print(plot)




results_table_fix <- results_table %>%
  mutate(ensgene = gsub("\\..*", "", Geneid)) %>%
  left_join(annotables::grch38[,c(1,3)]) %>%
  mutate(gene_name = symbol)
results_table_fix_2 <- results_table_fix[!duplicated(results_table_fix$gene_name),] #%>%
  #unique(symbol)

big_data_with_counts <- big_data_ce %>%
  left_join(results_table_fix_2) %>%
  dplyr::mutate(gene_color = ifelse(padj < 0.1 & log2FoldChange < 0, "down", 
                                    ifelse(padj < 0.1 & log2FoldChange > 0, "up", "not_sign"))) %>%
  dplyr::mutate(category = ifelse(Control_TDP43KD < 0.05 & Cycloheximide_TDP43KD > 0.1, "dark",
                                  ifelse(Control_TDP43KD > 0.1 & Cycloheximide_TDP43KD < 0.05, "anti-dark","else"))) 

#write.csv(big_data_with_counts, "antidark.csv")

big_data_with_counts2 <- big_data_with_counts[!is.na(big_data_with_counts$category),]
big_data_with_counts3 <- big_data_with_counts2[!is.na(big_data_with_counts2$gene_color),]
big_data_with_counts_long <- big_data_with_counts3 %>%
  pivot_longer(cols = c(Control_Control,
                        Cycloheximide_Control, 
                        Control_TDP43KD, 
                        Cycloheximide_TDP43KD))

big_data_with_counts_long$name <- factor(big_data_with_counts_long$name, 
                                levels = c("Control_Control",
                                           "Cycloheximide_Control", 
                                           "Control_TDP43KD", 
                                           "Cycloheximide_TDP43KD"))


plot <- big_data_with_counts_long %>%
  #dplyr::filter(value < 2 & value > -2) %>%
  #dplyr::filter(paste_into_igv_junction == "chr8:79611214-79616822") %>%
  ggplot(mapping = aes(x = name, y = value, group = paste_into_igv_junction, color = gene_color)) +
  geom_line(aes(), show.legend = T) +
  #geom_bar(fill ="#3288BD", stat = "identity", width = 0.5) +
  #geom_errorbar(aes(x = .id, ymin = mean_psi_per_lsv_junction - (stdev_psi_per_lsv_junction/sqrt(2)),
  #                           ymax = mean_psi_per_lsv_junction + (stdev_psi_per_lsv_junction/sqrt(2))), width = 0.3) +
  #geom_text_repel(aes(x = .id, y = mean_psi_per_lsv_junction, label = label_junction), point.padding = 0.3,
  #                nudge_y = 0.2, min.segment.length = 0.5, box.padding  = 2, max.overlaps = Inf, size=4, show.legend = F) +
  #scale_color_manual(values = c("gray","#3288BD",
  #                              "#D53E4F")) +
  #scale_alpha_manual(values = c(0.1,1)) +
  #xlab("TDP-43 protein levels") +
  facet_grid(facets = vars(category)) +
  ylab("PSI") +
  #scale_x_discrete(labels = c("100%", "77%", "25%", "8%", "4%", "0%")) +
  scale_y_continuous(breaks=seq(0,1,0.2),
                     labels = c("0%","20%","40%","60%","80%","100%")) +
  theme_classic()
print(plot)

big_data_with_counts3 %>%
  #dplyr::filter(gene_color == "down") %>%
  ggplot(aes(x = log2FoldChange)) +
  stat_ecdf(aes(color = gene_color)) +
  #geom_density(aes(fill = gene_color), alpha = 0.2) +
  #facet_grid(facets = vars(gene_color)) +
  theme_classic()




normed_counts_recap <- normed_counts %>%
  dplyr::mutate(Cycloheximide_Control = rowMeans(normed_counts[,c(2:5)])) %>%
  dplyr::mutate(Control_Control = rowMeans(normed_counts[,c(6:9)])) %>%
  dplyr::mutate(Cycloheximide_TDP43KD = rowMeans(normed_counts[,c(10:13)])) %>%
  dplyr::mutate(Control_TDP43KD = rowMeans(normed_counts[,c(14:17)])) %>%
  dplyr::mutate(ensgene = gsub("\\..*", "", Geneid)) %>%
  dplyr::select(18:22) %>%
  left_join(annotables::grch38[,c(1,3)])

results_table_fix %>%
  dplyr::group_by(symbol) %>%
  dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
  dplyr::filter(n > 1L) 

normed_counts_recap_long <- normed_counts_recap %>%
  pivot_longer(cols = c("Control_Control",
                        "Cycloheximide_Control", 
                        "Control_TDP43KD", 
                        "Cycloheximide_TDP43KD"))
normed_counts_recap_long$name <- factor(normed_counts_recap_long$name, 
                                levels = c("Control_Control",
                                           "Cycloheximide_Control", 
                                           "Control_TDP43KD", 
                                           "Cycloheximide_TDP43KD"))

  

plot <- normed_counts_recap_long %>%
  #dplyr::filter(value < 2 & value > -2) %>%
  #dplyr::filter(paste_into_igv_junction == "chr8:79611214-79616822") %>%
  ggplot(mapping = aes(x = name, y = value, group = symbol)) +
  geom_line(aes(), show.legend = F) +
  #geom_bar(fill ="#3288BD", stat = "identity", width = 0.5) +
  #geom_errorbar(aes(x = .id, ymin = mean_psi_per_lsv_junction - (stdev_psi_per_lsv_junction/sqrt(2)),
  #                           ymax = mean_psi_per_lsv_junction + (stdev_psi_per_lsv_junction/sqrt(2))), width = 0.3) +
  #geom_text_repel(aes(x = .id, y = mean_psi_per_lsv_junction, label = label_junction), point.padding = 0.3,
  #                nudge_y = 0.2, min.segment.length = 0.5, box.padding  = 2, max.overlaps = Inf, size=4, show.legend = F) +
  #scale_color_manual(values = c("gray","#3288BD",
  #                              "#D53E4F")) +
  #scale_alpha_manual(values = c(0.1,1)) +
  #xlab("TDP-43 protein levels") +
  ylab("PSI") +
  #scale_x_discrete(labels = c("100%", "77%", "25%", "8%", "4%", "0%")) +
  #scale_y_continuous(breaks=seq(0,1,0.2),
  #                   labels = c("0%","20%","40%","60%","80%","100%")) +
  theme_classic()
print(plot)


###add de to majiq output
#big_delta_narrow <- big_delta[,c(1,7,8,10,11,13,14,18,22,23,24,36)]

big_data_narrow <- big_data[,c(1,2,5,9,11,17)] %>%
  #group_by(lsv_junc) %>%
  #filter(n() == 4) %>%
  pivot_wider(names_from = .id, values_from = mean_psi_per_lsv_junction, values_fn = max) %>%
  filter(Control_Control < 0.05) %>%
  filter(Control_TDP43KD - Control_Control > 0.1 | Cycloheximide_TDP43KD - Control_Control > 0.1) %>%
  mutate(color_gene_name = as.character(ifelse(Cycloheximide_TDP43KD - Control_TDP43KD > 0.05, 
                                               "Delta PSI > 0.05", "Delta PSI < 0.05")))
names(big_data_narrow)[1] <- "symbol"
table(big_data_narrow$color_gene_name)

normed_counts_avg <- normed_counts_all[,c(2:18)] %>%
  left_join(annotables::grch38[,c(1,3)]) %>%
  mutate(Control_Control_de = rowMeans(.[,c(5:8)])) %>%
  mutate(Cycloheximide_Control_de = rowMeans(.[,c(1:4)])) %>%
  mutate(Control_TDP43KD_de = rowMeans(.[,c(13:16)])) %>%
  mutate(Cycloheximide_TDP43KD_de = rowMeans(.[,c(9:12)])) %>%
  select(c(18:22))

big_data_de_normed_counts <- left_join(big_data_narrow, normed_counts_avg) %>% 
  filter(!is.na(color_gene_name))
#write.csv(big_data_de, "~/Desktop/annotation_nmd.csv", row.names = F)

big_data_de_normed_counts %>%
  ggplot(aes(x = log2(Control_TDP43KD_de/Control_Control_de), 
             y = Control_TDP43KD-Control_Control,
             color = color_gene_name)) +
  geom_point() +
  geom_smooth(method = 'lm', se = F) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  #facet_wrap(facets = vars(color_gene_name)) +
  xlim(c(-4,4)) +
  ylim(c(0,1)) +
  scale_color_brewer(palette = "Set1") +
  #coord_fixed() +
  theme_classic()


big_data_de_results_table <- left_join(big_data_narrow, results_table_ctrl[,c(2,3,7,9)])
##trying to visualize CE-dependent DE 
big_data_de_results_table %>%
  filter(de_novo_junctions == 1 & Control_TDP43KD > 0.1 &
           padj < 0.05 & !is.na(color_gene_name)) %>%
  ggplot(aes(x = log2FoldChange, 
             y = Control_TDP43KD-Control_Control, color = color_gene_name)) +
  geom_point() +
  geom_smooth(method = 'lm', se = F) +
  facet_wrap(facets = vars(color_gene_name)) +
  theme_classic()


big_delta_narrow <- big_delta[,c(1,7,8,10,11,13,14,18,22,23,24,36)] %>%
  filter(.id == "ControlControl-ControlTDP43KD" & base_mean_psi < 0.05 & probability_changing > 0.9) %>%
  names(big_delta_narrow)[2] <- "symbol"

big_delta_de_results_table <- left_join(big_delta_narrow, results_table_ctrl[,c(2,3,7,9)])
big_delta_de_results_table %>%
  ggplot(aes(x = log2FoldChange, 
             y = mean_dpsi_per_lsv_junction)) +
  geom_point() +
  geom_smooth(method = 'lm', se = F) +
  #facet_wrap(facets = vars(color_gene_name)) +
  theme_classic()
