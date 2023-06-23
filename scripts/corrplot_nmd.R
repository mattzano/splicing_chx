##nmd prediction stats

## import nmd prediction table
non_nmd_expected <- read.table("data/non_nmd_expected.csv", sep = ",", header = T)
## compare to slope plots
big_data_rescue <- big_data_chx %>%
  group_by(lsv_junc) %>%
  dplyr::filter(n() == 4) %>% #???
  dplyr::filter(de_novo_junctions == 1) %>%
  dplyr::filter(mean_psi_per_lsv_junction[.id == input_list[1]] < 0.05) %>%
  dplyr::filter(mean_psi_per_lsv_junction[.id == input_list[2]] < 0.05) %>%
  dplyr::filter(mean_psi_per_lsv_junction[.id == input_list[3]] > 0.05) %>% 
  dplyr::filter(mean_psi_per_lsv_junction[.id == input_list[4]] > 0.1) %>%
  dplyr::filter((mean_psi_per_lsv_junction[.id == input_list[3]] - mean_psi_per_lsv_junction[.id == input_list[1]] > 0.05 |
                   mean_psi_per_lsv_junction[.id == input_list[4]] - mean_psi_per_lsv_junction[.id == input_list[1]] > 0.05)) %>%
  #dplyr::filter(mean_psi_per_lsv_junction[.id == input_list[4]] - mean_psi_per_lsv_junction[.id == input_list[3]] > 0) %>%# &
  mutate(color_gene_name = as.character(ifelse(mean_psi_per_lsv_junction[.id == input_list[4]] - mean_psi_per_lsv_junction[.id == input_list[3]] > 0.05, 
                                               "NMD rescued", "non-NMD rescued"))) %>%
  dplyr::filter(.id == input_list[3])
  #mutate(alpha_gene_name = as.character(ifelse(mean_psi_per_lsv_junction[.id == input_list[4]] > 0.2 &
  #                                               gene_name %in% gene_list, 1, 0))) %>%#c("chr19:17641556-17642414","chr8:79611214-79616822"),"chr20:63439708-63444659","chr10:3099819-3101365","chr2:241668985-241670726","chr19:7169831-7170537"), 
  #mutate(label_junction = case_when(.id == input_list[4] & mean_psi_per_lsv_junction[.id == input_list[4]] > 0.2 &
  #                                    gene_name %in% gene_list ~ gene_name, T ~ ""))

non_nmd_table <- right_join(non_nmd_expected, big_data_rescue, by = "paste_into_igv_junction") %>%
  distinct(paste_into_igv_junction, .keep_all = T) %>%
  mutate(PTC_corrected = ifelse(is.na(PTC), "PTC predicted", "non-PTC predicted"))

table(non_nmd_table$PTC_corrected, non_nmd_table$color_gene_name)

non_nmd_table %>%
  ggplot(aes(x = color_gene_name, y = mean_psi_per_lsv_junction, color = PTC_corrected)) +
  geom_point(position = position_dodge2(width = 0.8)) +
  theme_classic()

###try and get cleaner data from manual curation
manual_validation_sy5y <- read_xlsx("~/Documents/phd/research_lines/tdp-43 curves/manual_validation_curves/validation_curves.xlsx") %>%
  mutate(paste_into_igv_junction = paste0(chr, ":", start, "-", end))

manual_validated_nmd <- left_join(non_nmd_table, manual_validation_sy5y, by = "paste_into_igv_junction") %>%
  dplyr::filter(!is.na(Type))

manual_validated_nmd %>%
  ggplot(aes(x = color_gene_name, y = mean_psi_per_lsv_junction, color = PTC_corrected)) +
  geom_point(position = position_dodge2(width = 0.8)) +
  theme_classic()

library(corrplot)
chisq <- chisq.test(table(manual_validated_nmd$PTC_corrected, manual_validated_nmd$color_gene_name))
M <- as.matrix(chisq$residuals)
png(file = "~/Desktop/corr_nmd_ptc.png")
corrplot(M, is.corr = F, cl.pos = 'n')
dev.off()

