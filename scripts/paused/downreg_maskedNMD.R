
#names(normed_counts)[1] <- "ensgene"
#normed <- normed_counts %>%
#  mutate(ensgene = gsub("\\..*", "", Geneid)) %>%
#downreg <- normed %>%  
#left_join(annotables::grch38[,c(1,3)], by = "ensgene")

#downreg <- normed %>%
#  mutate(ctrl = rowMeans(dplyr::select(normed, contains("CTRL_ctrl")), na.rm = TRUE)) %>%
#  #mutate(stDev = apply(.[(2:5)],1,sd)) ##
#  mutate(ctrl_sd = apply(.[(2:5)],1,sd)) %>%
#  mutate(tdp = rowMeans(dplyr::select(normed, contains("DOX_ctrl")), na.rm = TRUE)) %>%
#  mutate(tdp_sd = apply(.[(6:9)],1,sd)) %>%
#  mutate(delta = tdp / ctrl) %>%
#  mutate(delta_log2 = log2(delta)) %>%
#  mutate(t_test = (abs(tdp - ctrl)) / (sqrt((tdp_sd^2 / 4) + (ctrl_sd^2 / 4)))) %>% ### welch
#  mutate(df = ((tdp_sd^2 / 4) + (ctrl_sd^2 / 4))^2 / ((((tdp_sd^2 / 4)^2) /3) + (((ctrl_sd^2 / 4)^2) /3)) ) %>%
#  mutate(p = pt(t_test, df, lower.tail = F)) %>%
# mutate(color = ifelse(abs(delta_log2) > 1 & p < 0.01, 1, 0))

downreg <- results_table %>%
  mutate(ensgene = gsub("\\..*", "", Geneid)) %>%
  left_join(annotables::grch38[,c(1,3)], by = "ensgene") %>%
  mutate(color = ifelse(abs(log2FoldChange) > 1 & padj < 0.01, 1, 0))
  
write.csv(downreg, "downreg.csv")


downreg %>%
  ggplot(aes(y = -log10(padj), x = log2FoldChange, color = color)) +
  geom_point()
  
downreg_ok <- downreg %>%
  dplyr::filter(log2FoldChange < 0 & color == 1)
names(downreg_ok)[9] <- "gene_name"

dark_cryptome <- big_data_filtereds %>%
  group_by(lsv_junc) %>%
  dplyr::filter(n() == 4) %>%
  dplyr::filter(mean_psi_per_lsv_junction[.id == "Cycloheximide_TDP43KD"] - mean_psi_per_lsv_junction[.id == "Control_Control"] > 0.05)

write.csv(dark_cryptome, "dark_crpytome.csv")


dark_cryp_down <- dark_cryptome %>%
  dplyr::filter(.id == "Cycloheximide_TDP43KD") %>%
  inner_join(downreg_ok, by = "gene_name")

dark_cryp_down %>%
  ggplot(aes(x = log2FoldChange, y = mean_psi_per_lsv_junction)) +
  geom_text_repel(aes(label = gene_name), data = dark_cryp_down[dark_cryp_down$mean_psi_per_lsv_junction > 0.1, ]) +
  geom_point()

write.csv(dark_cryp_down, "dark_cryptome_downreg.csv")


downreg_long <- downreg %>%
  pivot_longer(cols = c("ctrl","ctrl_sd", "tdp", "tdp_sd"), 
             names_to = "", values_to = "")


count_table <- function(featureCounts) {
  avg_read_counts <- featureCounts %>%
    mutate(dox_0125 = rowMeans(dplyr::select(featureCounts, contains("DOX_0.0125")), na.rm = TRUE)) %>%
    mutate(dox_0187 = rowMeans(dplyr::select(featureCounts, contains("DOX_0.0187")), na.rm = TRUE)) %>%
    mutate(dox_021 = rowMeans(dplyr::select(featureCounts, contains("DOX_0.021")), na.rm = TRUE)) %>%
    mutate(dox_025 = rowMeans(dplyr::select(featureCounts, contains("DOX_0.025")), na.rm = TRUE)) %>%
    mutate(dox_075 = rowMeans(dplyr::select(featureCounts, contains("DOX_0.075")), na.rm = TRUE)) %>%
    mutate(no_dox = rowMeans(dplyr::select(featureCounts, contains("NT_0")), na.rm = TRUE)) %>%
    dplyr::select(c(20, 26, 21:25)) %>%
    pivot_longer(cols = c("no_dox","dox_0125", "dox_0187", "dox_021", "dox_025", "dox_075"), 
                 names_to = ".id", values_to = "avg_count") %>%
    group_by(gene_name, .id) %>%
    summarise(avg_count = max(avg_count))
  avg_read_counts$.id <- factor(avg_read_counts$.id, levels = input_list)