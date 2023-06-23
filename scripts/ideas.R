big_data_filtereds_exp <- big_delta_chx %>%
  dplyr::filter(.id == "ControlTDP43KD-CycloheximideTDP43KD") %>%
  mutate(color_gene_name = as.character(ifelse(mean_dpsi_per_lsv_junction > 0.05, 
                                               "NMD rescued", "non-NMD rescued")))
big_data_filtereds_exp <- big_data_chx %>%
  group_by(lsv_junc) %>%
  dplyr::filter(n() == 4) %>% #???
  #dplyr::filter(de_novo_junctions == 1) %>%
  #dplyr::filter(mean_psi_per_lsv_junction[.id == input_list[1]] < 0.05) %>%
  #dplyr::filter(mean_psi_per_lsv_junction[.id == input_list[2]] < 0.05) %>%
  #dplyr::filter(mean_psi_per_lsv_junction[.id == input_list[3]] < 0.05) %>% 
  #dplyr::filter(mean_psi_per_lsv_junction[.id == input_list[4]] > 0.1) %>%
  dplyr::filter((mean_psi_per_lsv_junction[.id == input_list[3]] - mean_psi_per_lsv_junction[.id == input_list[1]] > 0.05 |
                   mean_psi_per_lsv_junction[.id == input_list[4]] - mean_psi_per_lsv_junction[.id == input_list[1]] > 0.05)) %>%
  #dplyr::filter(mean_psi_per_lsv_junction[.id == input_list[4]] - mean_psi_per_lsv_junction[.id == input_list[3]] > 0) %>%# &
  mutate(color_gene_name = as.character(ifelse(mean_psi_per_lsv_junction[.id == input_list[4]] - mean_psi_per_lsv_junction[.id == input_list[3]] > 0.05, 
                                               "NMD rescued", "non-NMD rescued")))
  
write.table(big_data_filtereds_exp, "~/Desktop/nmd_chx.csv", sep = ",", row.names = F)



#### what if comparing tdp43 binding sites with genomic and transcriptomics information
## ie conservation / cdc scores - genomic
## gene location - 5'utr - intron - exon - 3'utr - compare with ule's papers
## do these info change upon different cell lines??

#load iclip data/postar data

#load conservation score data

#load gene info - relative to genomic location

#make scatter/box plots



