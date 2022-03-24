slope_plot_nmd <- function(gene_list){
big_data_filtereds <- big_data %>%
  group_by(lsv_junc) %>%
  filter(n() == 4) %>% #???
  filter(mean_psi_per_lsv_junction[.id == "Control_Control"] < 0.05) %>%
  filter(mean_psi_per_lsv_junction[.id == "Control_TDP43KD"] > 0.1) %>%
  #filter(mean_psi_per_lsv_junction[.id == "Cycloheximide_TDP43KD"] - mean_psi_per_lsv_junction[.id == "Control_TDP43KD"] > 0) %>% #&  #these tow filters are to make plot prettier
         #mean_psi_per_lsv_junction[.id == "Cycloheximide_TDP43KD"] - mean_psi_per_lsv_junction[.id == "Control_TDP43KD"] < 0.1) %>%
  filter((mean_psi_per_lsv_junction[.id == "Control_TDP43KD"] - mean_psi_per_lsv_junction[.id == "Control_Control"] > 0.05 |
           mean_psi_per_lsv_junction[.id == "Cycloheximide_TDP43KD"] - mean_psi_per_lsv_junction[.id == "Control_Control"] > 0.05)) %>% # &
           #mean_psi_per_lsv_junction[.id == "Cycloheximide_TDP43KD"] - mean_psi_per_lsv_junction[.id == "Cycloheximide_Control"] > 0.05) %>%   ###???
  mutate(color_gene_name = as.character(ifelse(mean_psi_per_lsv_junction[.id == "Cycloheximide_TDP43KD"] - mean_psi_per_lsv_junction[.id == "Control_TDP43KD"] > 0.05, 
                                               "Delta PSI > 0.05", "Delta PSI < 0.05"))) %>%
  mutate(alpha_gene_name = as.character(ifelse(gene_name %in% gene_list, 1, 0))) %>%#c("chr19:17641556-17642414","chr8:79611214-79616822"),"chr20:63439708-63444659","chr10:3099819-3101365","chr2:241668985-241670726","chr19:7169831-7170537"), 
  mutate(label_junction = case_when(.id =="Control_TDP43KD" & gene_name %in% gene_list ~ gene_name, T ~ ""))

big_data_filtereds %>%
  filter(.id %in% c("Control_TDP43KD","Cycloheximide_TDP43KD")) %>%
  ggplot(mapping = aes(x = .id, y = mean_psi_per_lsv_junction)) +
  facet_wrap(facets = vars(color_gene_name)) +
  geom_point(aes(color = color_gene_name, alpha = alpha_gene_name, group = lsv_junc), show.legend = F) + 
  geom_line(aes(color = color_gene_name, alpha = alpha_gene_name, group = lsv_junc), show.legend = F) +
  geom_text_repel(aes(x = .id, y = mean_psi_per_lsv_junction, label = label_junction), point.padding = 0.3,
                  nudge_y = 0.2, min.segment.length = 0.5, box.padding  = 2, max.overlaps = Inf, size=4, show.legend = F) +
  scale_color_manual(values = c(#"#ABDDA4", 
                                "#3288BD",
                                "#D53E4F")) +
  scale_alpha_manual(values = c(0.02,1)) +
  xlab("") +
  ylab("PSI") +
  scale_x_discrete(labels = c("TDP43-KD", "TDP43-KD+CHX")) +
  theme_classic()
}

#avg_read_counts <- featureCounts %>%
#  mutate(Control_Control = rowMeans(dplyr::select(featureCounts, contains("CTRL_ctrl")), na.rm = TRUE)) %>%
#  mutate(Cycloheximide_Control = rowMeans(dplyr::select(featureCounts, contains("CTRL_chx")), na.rm = TRUE)) %>%
#  mutate(Control_TDP43KD = rowMeans(dplyr::select(featureCounts, contains("DOX_ctrl")), na.rm = TRUE)) %>%
#  mutate(Cycloheximide_TDP43KD = rowMeans(dplyr::select(featureCounts, contains("DOX_chx")), na.rm = TRUE)) %>%
#  dplyr::select(c(18:22)) %>%
#  pivot_longer(cols = c("Control_Control", "Cycloheximide_Control", "Control_TDP43KD", "Cycloheximide_TDP43KD"), 
#               names_to = ".id", values_to = "avg_count") %>%
#  group_by(gene_name, .id) %>%
#  summarise(avg_count = max(avg_count))

#big_data_avg_counts <- left_join(big_data, avg_read_counts)
#big_data_avg_counts_annot <- left_join(big_data_avg_counts, big_delta_filter)

#big_data_avg_counts_annot %>%
#  filter(.id %in% c("Control_TDP43KD","Cycloheximide_TDP43KD")) %>%
#  ggplot(mapping = aes(x = .id, y = avg_count)) +
  #facet_wrap(facets = vars(junc_cat)) +
#  geom_point(aes(color = color_gene_name, alpha = alpha_gene_name, group = lsv_junc), show.legend = T) + 
  #geom_violin(aes(x = .id, y = mean_psi_per_lsv_junction_normal_fake, group = color_gene_name, color = color_gene_name), show.legend = F) +
#  geom_line(aes(color = color_gene_name, alpha = alpha_gene_name, group = lsv_junc), show.legend = F) +
#  geom_text_repel(aes(label = label_junction), point.padding = 0.3,
#                  nudge_y = 0.2, min.segment.length = 0.5, box.padding  = 2, max.overlaps = Inf, size=4, show.legend = F) +
#  scale_color_manual(values = c(#"#ABDDA4", 
#                                "#3288BD",
#                                "#D53E4F")) +
#  scale_alpha_manual(values = c(0.1,1)) +
  #scale_x_log10() +
#  scale_y_log10() +
#  theme_classic()