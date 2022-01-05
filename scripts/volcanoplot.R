if (!require("pacman")) install.packages("pacman")
library(pacman)
p_load("data.table",
       "ggnewscale",
       "tidyverse",
       "ggplot2",
       "rcompanion",
       "lme4",
       "ggpubr",
       "forcats",
       "scales",
       "here",
       "ggthemes",
       "ggrepel",
       "DESeq2",
       "patchwork",
       "rstatix",
       "ggeasy")

ipsc_splicing <- fread(file.path(here::here(),"chx_splicing_results.csv"))
ipsc_splicing <- separate_rows(ipsc_splicing, c(mean_dpsi_per_lsv_junction, probability_changing, probability_non_changing, ControlControl_mean_psi, ControlTDP43KD_mean_psi, lsv_type, de_novo_junctions, junctions_coords), 
                               sep = c(";"), convert = T)

splicing_dots_tables <- ipsc_splicing %>%
  mutate(junction_name = case_when(gene_name %in% c("UNC13A","AGRN",
                                                    "UNC13B","PFKP","SETD5",
                                                    "ATG4B","STMN2") &
                                     probability_changing > 0.9  &
                                    mean_dpsi_per_lsv_junction > 0 ~ gene_name,
                                   T ~ "")) %>%
  mutate(`Novel Junction` = de_novo_junctions == 0) %>%
  mutate(log10_test_stat = -log10(1 - probability_changing)) %>%
  mutate(log10_test_stat = ifelse(is.infinite(log10_test_stat), 4.5, log10_test_stat)) %>%
  mutate(graph_alpha = ifelse(probability_changing > 0.9, 1, 0.2)) %>%
  mutate(label_junction = case_when(gene_name %in% c("UNC13A","AGRN",
                                                     "UNC13B","PFKP","SETD5",
                                                     "ATG4B","STMN2") &
                                      probability_changing > 0.9  &
                                     mean_dpsi_per_lsv_junction > 0 ~ junction_name,
                                    T ~ ""))


fig1a <- ggplot() +
  geom_point(data = splicing_dots_tables %>% filter(de_novo_junctions != 0),
             aes(x = mean_dpsi_per_lsv_junction, y = log10_test_stat,
                 alpha = graph_alpha, fill = "Annotated Junction"), pch = 21, size = 2) +
  geom_point(data = splicing_dots_tables %>% filter(de_novo_junctions == 0),
             aes(x = mean_dpsi_per_lsv_junction, y = log10_test_stat,
                 alpha = graph_alpha, fill = "Novel Junction"), pch = 21, size = 2) +
  geom_text_repel(data = splicing_dots_tables[19 != ""],
                  aes(x = mean_dpsi_per_lsv_junction, y =log10_test_stat,
                      label = label_junction,
                      color = as.character(de_novo_junctions)), point.padding = 0.3,
                  nudge_y = 0.2,
                  min.segment.length = 0,
                  box.padding  = 2,
                  max.overlaps = Inf,
                  size=4, show.legend = F) +
  geom_hline(yintercept = -log10(1 - .9)) +
  geom_vline(xintercept = -0,linetype="dotted") +
  scale_fill_manual(name = "",
                    breaks = c("Annotated Junction", "Novel Junction"),
                    values = c("Annotated Junction" = "#648FFF", "Novel Junction" = "#fe6101") ) +
  scale_color_manual(name = "",
                     breaks = c("0", "1"),
                     values = c("1" = "#648FFF", "0" = "#fe6101") ) +
  guides(alpha = FALSE, size = FALSE) +
  theme(legend.position = 'top') +
  ggpubr::theme_pubr() +
  xlab("delta PSI") +
  ylab(expression(paste("-Lo", g[10], " Test Statistic"))) +
  theme(text = element_text(size = 24)) +
  theme(legend.text=element_text(size=22)) +
  xlim(-1,1) +
  scale_x_continuous(labels = scales::percent)
fig1a
