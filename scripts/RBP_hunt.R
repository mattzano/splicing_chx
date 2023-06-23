###cryptics in chx alone


####rbp with CEs --- either using scoring, Yeo 2020 paper, or genes with RRM from uniprot!!

chx_cryptics <- big_delta %>%
  dplyr::filter(.id == "ControlControl-CycloheximideControl" | .id == "ControlControl-ControlTDP43KD") %>%
  group_by(lsv_junc) %>%
  dplyr::filter(mean_dpsi_per_lsv_junction > 0.1 &
                  base_mean_psi < 0.05) %>%
  dplyr::filter(n() == 2)

#chx_cryptics <- big_data %>%
#  dplyr::filter(.id != "Cycloheximide_TDP43KD") %>%
#  group_by(lsv_junc) %>%
#  dplyr::filter(n() == 3) %>%
#  dplyr::filter(mean_psi_per_lsv_junction[.id == "Control_Control"] < 0.05 &
#                  mean_psi_per_lsv_junction[.id == "Cycloheximide_Control"] > 0.1 &
#                  mean_psi_per_lsv_junction[.id == "Control_TDP43KD"] > 0.1) %>%
#  dplyr::filter(n() == 3)

#table(chx_cryptics$junc_cat)
#table(big_delta$.id)

rbp_list <- read.table("~/Desktop/rbp_bed/Table_HS_RBP.txt", header = T)
#####als_gene <- read.csv("~/Desktop/Supplementary File 1.csv")

tdp_cryptics <- big_data %>%
  dplyr::filter(.id != "Cycloheximide_TDP43KD" & .id != "Cycloheximide_Control") %>%
  group_by(lsv_junc) %>%
  dplyr::filter(n() == 2) %>%
  dplyr::filter(mean_psi_per_lsv_junction[.id == "Control_Control"] < 0.05 &
                  #mean_psi_per_lsv_junction[.id == "Cycloheximide_Control"] > 0.1 &
                  mean_psi_per_lsv_junction[.id == "Control_TDP43KD"] > 0.1) %>%
  dplyr::filter(n() == 2)

tdp_cryptics <- big_delta %>%
  dplyr::filter(.id == "ControlControl-ControlTDP43KD") %>%
  group_by(lsv_junc) %>%
  dplyr::filter(mean_dpsi_per_lsv_junction > 0.1 &
                  base_mean_psi < 0.05)


tdp_rbp <- tdp_cryptics %>%
  dplyr::filter(gene_name %in% rbp_list$Gene_Name &
                  probability_changing > 0.9)#[rbp_list$RBP2GO_Score > 20])

als_cryptics <- tdp_cryptics %>%
  dplyr::filter(gene_name %in% als_gene$gene_name)


patient_junction_list <- patient_junction_list[,2]

patient_junction_list <- patient_junction_list %>%
  pull()

big_data_patient <- big_data %>%
  dplyr::filter(paste_into_igv_junction %in% patient_junction_list)

big_delta_patient <- big_delta %>%
  dplyr::filter(paste_into_igv_junction %in% patient_junction_list)


big_delta_patient_filtered <- big_delta_patient %>%
  dplyr::filter(.id == "ControlControl-ControlTDP43KD" | .id == "ControlControl-CycloheximideTDP43KD") %>%
  dplyr::filter(base_mean_psi < 0.05)


brain_junction_list <- big_delta %>%
  dplyr::filter(.id == "ControlControl-ControlTDP43KD" | .id == "ControlControl-CycloheximideTDP43KD") %>%
  dplyr::filter(base_mean_psi < 0.05 & mean_dpsi_per_lsv_junction > 0.1) %>%
  pull(paste_into_igv_junction)
