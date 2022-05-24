big_data_conc <- read_csv("~/Desktop/big_data_conc.csv", 
                          col_types = cols(.id = col_factor(levels = c("no_dox", 
                                                            "dox_0125", "dox_0187", "dox_021", 
                                                            "dox_025", "dox_075"))))
big_delta_conc <- read.csv("~/Desktop/big_delta_conc.csv")
big_data_chx <- big_data
big_delta_chx <- big_delta

big_delta_chx <- big_delta_chx %>%
  dplyr::filter(junc_cat == "novel_acceptor" | junc_cat == "novel_combo" | junc_cat == "novel_donor")
big_delta_conc <- big_delta_conc %>%
  dplyr::filter(junc_cat == "novel_acceptor" | junc_cat == "novel_combo" | junc_cat == "novel_donor")


big_data_chx_filtereds <- big_data_chx %>%
  dplyr::filter(lsv_junc %in% big_delta_chx$lsv_junc) %>%
  group_by(lsv_junc) %>%
  dplyr::filter(n() == 4) %>%
  dplyr::filter(mean_psi_per_lsv_junction[.id == "Control_Control"] < 0.05) %>%
  dplyr::filter(mean_psi_per_lsv_junction[.id == "Control_TDP43KD"] > 0.1) %>%
  dplyr::filter((mean_psi_per_lsv_junction[.id == "Control_TDP43KD"] - mean_psi_per_lsv_junction[.id == "Control_Control"] > 0.05 |
                 mean_psi_per_lsv_junction[.id == "Cycloheximide_TDP43KD"] - mean_psi_per_lsv_junction[.id == "Control_Control"] > 0.05)) %>% # &
  dplyr::mutate(color_gene_name = as.character(ifelse(mean_psi_per_lsv_junction[.id == "Cycloheximide_TDP43KD"] - mean_psi_per_lsv_junction[.id == "Control_TDP43KD"] > 0.05, 
                                            "nmd", "not_nmd"))) %>%
  dplyr::mutate(mean_dpsi = mean_psi_per_lsv_junction[.id == "Cycloheximide_TDP43KD"] - mean_psi_per_lsv_junction[.id == "Control_TDP43KD"])
  #mutate(alpha_gene_name = as.character(ifelse(gene_name %in% gene_list, 1, 0))) %>%#c("chr19:17641556-17642414","chr8:79611214-79616822"),"chr20:63439708-63444659","chr10:3099819-3101365","chr2:241668985-241670726","chr19:7169831-7170537"), 
  #mutate(label_junction = case_when(.id =="Control_TDP43KD" & gene_name %in% gene_list ~ gene_name, T ~ ""))
  
#big_data_fused <- big_data_conc %>%
#  dplyr::mutate(chx = as.factor(as.character(ifelse(lsv_junc %in% big_data_chx$lsv_junc[big_data_chx_filtereds$color_gene_name == "nmd"], "nmd", 
#                                                    ifelse(lsv_junc %in% big_data_chx$lsv_junc[big_data_chx_filtereds$color_gene_name == "not_nmd"], "not_nmd", "none")))))
  
big_data_conc_filtereds <- big_data_conc %>%
  dplyr::mutate(chx = as.factor(as.character(ifelse(lsv_junc %in% big_data_chx_filtereds$lsv_junc[big_data_chx_filtereds$color_gene_name == "nmd"], "nmd", 
                                                    ifelse(lsv_junc %in% big_data_chx_filtereds$lsv_junc[big_data_chx_filtereds$color_gene_name == "not_nmd"], "not_nmd", "none"))))) %>%
  dplyr::filter(lsv_junc %in% big_delta_conc$lsv_junc) %>% #try without cryptics
  group_by(lsv_junc) %>% 
  dplyr::filter(n() == 6) %>%
  dplyr::filter(mean_psi_per_lsv_junction[.id == "no_dox"]  < 0.05) %>% 
  dplyr::filter(mean_psi_per_lsv_junction[.id == "dox_075"] > 0.1) %>%
  dplyr::filter((mean_psi_per_lsv_junction[.id == "dox_075"]  - mean_psi_per_lsv_junction[.id == "dox_025"]  > - 0.05) &
                (mean_psi_per_lsv_junction[.id == "dox_025"]  - mean_psi_per_lsv_junction[.id == "dox_021"]  > 0) &
                (mean_psi_per_lsv_junction[.id == "dox_021"]  - mean_psi_per_lsv_junction[.id == "dox_0187"] > 0) &
                (mean_psi_per_lsv_junction[.id == "dox_0187"] - mean_psi_per_lsv_junction[.id == "dox_0125"] > 0) &
                (mean_psi_per_lsv_junction[.id == "dox_0125"] - mean_psi_per_lsv_junction[.id == "no_dox"]   > - 0.05)) %>%
  dplyr::mutate(color_gene_name_fused = as.factor(as.character(ifelse(mean_psi_per_lsv_junction[.id == "dox_075"] > 0.1 & mean_psi_per_lsv_junction[.id == "dox_025"] < 0.1, "late", 
                                                                ifelse(mean_psi_per_lsv_junction[.id == "dox_0125"] > 0.1, "early", "none"))))) %>%
  dplyr::mutate(alpha_gene_name = as.factor(as.character(ifelse(color_gene_name_fused == "early" | color_gene_name_fused == "late", 2, 1)))) %>%
  dplyr::mutate(label_junction = case_when((.id =="dox_0187" & color_gene_name_fused == "early") | (.id == "dox_075" & color_gene_name_fused == "late") ~ gene_name, T ~ ""))
  
  
big_data_conc_filtereds %>%
  dplyr::filter(chx == "nmd" | chx == "not_nmd") %>%
  ggplot(mapping = aes(x = .id, y = mean_psi_per_lsv_junction, group = lsv_junc)) +
  geom_line(aes(color = chx, alpha = alpha_gene_name), show.legend = T) +
  geom_hline(aes(yintercept = 0.1)) +
  geom_text_repel(aes(x = .id, y = mean_psi_per_lsv_junction, label = label_junction), point.padding = 0.3,
                  nudge_y = 0.2, min.segment.length = 0.5, box.padding = 2, max.overlaps = Inf, size = 4, show.legend = F) +
  scale_color_manual(values = c(#"gray",
                                "#3288BD",
                                "#D53E4F")) +
  scale_alpha_manual(values = c(0.1,1)) +
  xlab("TDP-43 protein levels") +
  ylab("PSI") +
  scale_x_discrete(labels = c("100%","77%", "25%", "8%", "4%", "0%")) +
  scale_y_continuous(breaks=seq(0,1,0.2)) + #,
    #labels = c("0%","20%","40%","60%","80%","100%")) +
  theme_classic()
#print(plot)
