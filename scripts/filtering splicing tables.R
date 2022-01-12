gene_list <- c("ARHGAP32", "ITGA7", "CEP290", "PRELID3A", "DLGAP1", "UNC13A", "CELF5", "HDGFL2", "ATG4B", "KCNQ2",
               "KALRN", "CAMK2B", "ADCY1", "STMN2", "PHF2", "IGSF21", "TFAP2E", "ADARB2", "PFKP", "SYT7", "RSF1", 
               "G2E3", "KIAA0753", "CBARP", "INSR", "PXDN", "TRAPPC12", "SETD5", "TMEM175", "EPB41L4A", "MTRR", 
               "ACTL6B", "WASL", "IQCE", "TRRAP", "SEMA4D", "C1orf194", "ELAVL3", "ONECUT1", "KNDC1", "RBFOX3",
               "TLX1", "MAST1", "ZNF65", "PDE9A", "DAPK1", "COPS9", "AARS", "PSD", "PREX2")

input_splicing_ctrl_ctrltdp <- fread(file.path(here::here(), "data", "ControlControl-ControlTDP43KD.csv"))
input_splicing_ctrl_ctrltdp <- separate_rows(input_splicing_ctrl_ctrltdp, 
                                             c(names(input_splicing_ctrl_ctrltdp[,c(4:9,12,15)])), 
                                             sep  = c(";"), convert = T)
filtered_input_splicing_ctrl_ctrltdp <- input_splicing_ctrl_ctrltdp %>%
  filter(probability_changing > 0.9  &
           mean_dpsi_per_lsv_junction > 0 &
           de_novo_junctions == 1)

input_splicing_ctrl_ctrlchx <- fread(file.path(here::here(), "data", "ControlControl-CycloheximideControl.csv"))
input_splicing_ctrl_ctrlchx <- separate_rows(input_splicing_ctrl_ctrlchx, 
                                             c(names(input_splicing_ctrl_ctrlchx[,c(4:9,12,15)])), 
                                             sep  = c(";"), convert = T)
filtered_input_splicing_ctrl_ctrlchx <- input_splicing_ctrl_ctrlchx %>%
  filter(gene_name %in% gene_list &
    probability_changing > 0.9  &
           mean_dpsi_per_lsv_junction > 0)

input_splicing_tdp_ctrlchx <- fread(file.path(here::here(), "data", "ControlTDP43KD-CycloheximideTDP43KD.csv"))
input_splicing_tdp_ctrlchx <- separate_rows(input_splicing_tdp_ctrlchx, 
                                             c(names(input_splicing_tdp_ctrlchx[,c(4:9,12,15)])), 
                                             sep  = c(";"), convert = T)
filtered_input_splicing_tdp_ctrlchx <- input_splicing_tdp_ctrlchx %>%
  filter(gene_name %in% gene_list &
    probability_changing > 0.9  &
           mean_dpsi_per_lsv_junction > 0)

input_splicing_chx_ctrltdp <- fread(file.path(here::here(), "data", "CycloheximideControl-CycloheximideTDP43KD.csv"))
input_splicing_chx_ctrltdp <- separate_rows(input_splicing_chx_ctrltdp, 
                                             c(names(input_splicing_chx_ctrltdp[,c(4:9,12,15)])), 
                                             sep  = c(";"), convert = T)
filtered_input_splicing_chx_ctrltdp <- input_splicing_chx_ctrltdp %>%
  filter(probability_changing > 0.9  &
           mean_dpsi_per_lsv_junction > 0 &
           de_novo_junctions == 1)


diff_filtered_input_splicing_chx_ctrltdp <- filtered_input_splicing_chx_ctrltdp %>%
  filter(!lsv_id %in% filtered_input_splicing_ctrl_ctrltdp$lsv_id)

diff_filtered_input_splicing_ctrl_ctrltdp <- filtered_input_splicing_ctrl_ctrltdp %>%
  filter(!lsv_id %in% filtered_input_splicing_chx_ctrltdp$lsv_id)
