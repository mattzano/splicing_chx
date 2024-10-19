#### plot differential expression for top 50 up and 50 down in upf1 kd - 
### then check the same genes in tdp43 and tdp43/upf1kd

## do the same between chx and upf1

##### first make the differential expression files for upf1 vs control
outputdir = "~/Documents/Github/splicing_chx/results/deseq2/"
metadata_dir = "~/Documents/Github/splicing_chx/data/metadata_upf1.csv"  
###make this
tx2gene = "~/Desktop/rbp_bed/gencode.v40.tx2gene.csv"
column_name = "treatment"
baseline = "ctrl"
contrast = "UPF1KD"
controls_name = "Control"
contrast_name = "UPF1KD"
#column_name = opt$column_name

output_path=paste0(outputdir,controls_name,"-",contrast_name,".")

# ============================ section 1: import data ===========================

tx2gene <- data.table::fread(tx2gene,header=FALSE)
tx2gene <- tx2gene[,c(1:2)]
colnames(tx2gene) = c("TXNAME", "GENEID")

#(1) First read in the metadata. if only a subset of the files are used, the opt$pattern option will be taken.

#salmon_deseq2 <- function(mutation, sample_dir, metadata_dir) {
metadata_dir = "~/Documents/Github/splicing_chx/data/metadata_upf1.csv"  
metadata = read.csv(metadata_dir) %>% 
    #dplyr::select(sample, !!as.symbol(column_name)) %>% 
    dplyr::filter(condition == 0)
sample_dir = "~/Documents/GitHub/splicing_chx/data/salmon_upf1/"
  
  #(2) Generate a vector of the wanted file names.
  files = unique(file.path(sample_dir,metadata$sample,"quant.sf")) 
  names(files) = unique(metadata$sample)
  
  
  #(3) To check if all the files exist
  if(all(file.exists(files)) == FALSE) {
    stop("It seems that I cannot find those files...Please check if your directory is correct.")
  }
  
  
  # ====================== section 3: import salmon files ==============================
  # files is a vector of directory where quant.sf file locates.
  # just ignore the version... to make it easier for following steps.
  txi.tx <- tximport(files, 
                     type="salmon", 
                     tx2gene=tx2gene,
                     ignoreTxVersion = TRUE,
                     ignoreAfterBar = TRUE,
                     txOut = TRUE)
  
  txi.sum <- summarizeToGene(txi.tx, tx2gene)
  
  keep <- rowSums(edgeR::cpm(txi.sum$counts) > 5) >= 2
  print("Filtered Genes by CPM greater than 5 in a least 2 samples")
  print(table(keep))
  txi.sum$counts <- txi.sum$counts[keep, ]
  txi.sum$abundance <- txi.sum$abundance[keep, ]
  txi.sum$length <- txi.sum$length[keep, ]
  
  # make it csv
  TPM_transcripts = as.data.frame(txi.tx$abundance) %>% 
    tibble::rownames_to_column(.,var="transcript_id")
  TPM_gene = as.data.frame(txi.sum$abundance) %>% 
    tibble::rownames_to_column(.,var="gene_id")
  
  #write.csv(TPM_transcripts,"TPM_transcripts.csv")
  #write.csv(TPM_gene,"TPM_gene.csv")
  
  
  # ========================================== section 4: RUN A DEFAULT DESEQ 2 (optional) =============================================================
  
  
  dds = DESeqDataSetFromTximport(txi.sum,
                                 colData = metadata,
                                 design = ~ treatment) 
  
  
  
  # 'Note that the tximport-to-DESeq2 approach uses estimated gene counts from the transcript abundance quantifiers, but not normalized counts' -- <Deseq2 vignette> (just a note - Deseq() wraps the normalization step inside)
  # perform the Deseq function
  dds = DESeq(dds)
  
  # Now, extract the result and named them by their contrast group
  results_table_upf1vsctrl <<- results(dds) %>%
    as.data.frame() %>%
    tibble::rownames_to_column('Geneid') %>%
    mutate(ensgene = gsub("\\..*", "", Geneid)) %>%
    left_join(annotables::grch38[c(1,3)])
  
  # Now, extract the DESeq2 normed counts
  normed_counts_upf1vsctrl <<- counts(dds, normalized = TRUE) %>%
    as.data.frame() %>%
    tibble::rownames_to_column('Geneid') %>%
    mutate(ensgene = gsub("\\..*", "", Geneid)) %>%
    left_join(annotables::grch38[c(1,3)])
  summary(normed_counts)
  #return(results_table)
  
write.table(results_table_upf1vsctrl, "~/Desktop/results_table_upf1_vs_ctrl.csv",
            sep = ",", quote = F, row.names = F)

results_table_upf1vsctrl_filtered <- results_table_upf1vsctrl %>%
  filter(log2FoldChange > 1 | log2FoldChange < -1) %>%
  filter(!is.na(symbol))

results_table_UPF1_ctrl_filtered <- results_table_UPF1_ctrl %>%
  filter(symbol %in% results_table_upf1vsctrl_filtered$symbol) %>%
  filter(!is.na(symbol))

results_table_upf1vsctrl_filtered_2 <- results_table_upf1vsctrl_filtered %>%
  filter(symbol %in% results_table_UPF1_ctrl_filtered$symbol)

ggplot(data = results_table_upf1vsctrl_filtered_2, aes(x = reorder(symbol, log2FoldChange), y = 0, fill = log2FoldChange)) +
  geom_tile() +
  geom_tile(data = results_table_UPF1_ctrl_filtered, aes(x = symbol, y = -1, fill = log2FoldChange)) +
  theme_classic() +
  colorspace::scale_fill_continuous_diverging() +
  theme(axis.text.x = element_text(angle = 90))
 
results_table_upf1vsctrl_filtered_2 %>% pull(symbol) %>% write.csv("~/Desktop/top50_upf1.csv")


###concordance differential expression
results_table_upf1_filtered2 <- results_table_UPF1_ctrl %>%
  dplyr::filter(log2FoldChange > 1.5 | log2FoldChange < -1.5) %>%
  dplyr::filter(!is.na(symbol))
results_table_chx_filtered <- results_table_cont %>%
  dplyr::filter(symbol %in% results_table_upf1_filtered2$symbol) %>%
  #dplyr::filter(log2FoldChange > 1 | log2FoldChange < -1) %>%
  dplyr::filter(!is.na(symbol))
results_table_upf1_filtered <- results_table_upf1_filtered2 %>%
  dplyr::filter(symbol %in% results_table_chx_filtered$symbol)

ggplot(data = results_table_upf1_filtered, aes(x = reorder(symbol, log2FoldChange), y = 0, fill = log2FoldChange)) +
  geom_tile() +
  geom_tile(data = results_table_chx_filtered, aes(x = symbol, y = -1, fill = log2FoldChange)) +
  theme_classic() +
  colorspace::scale_fill_continuous_diverging() +
  theme(axis.text.x = element_text(angle = 90))

###concordance splicing
input_list <- c("Control_Control",
                "Cycloheximide_Control", 
                "Control_TDP43KD", 
                "Cycloheximide_TDP43KD")
big_data_chx_filtered <- big_data_chx %>%
  group_by(lsv_junc) %>%
  dplyr::filter(n() == 4) %>% #???
  dplyr::filter(de_novo_junctions == 1) %>%
  dplyr::filter(mean_psi_per_lsv_junction[.id == input_list[1]] < 0.05) %>%
  dplyr::filter(mean_psi_per_lsv_junction[.id == input_list[2]] < 0.05) %>%
  #dplyr::filter(mean_psi_per_lsv_junction[.id == input_list[3]] < 0.05) %>% 
  dplyr::filter(mean_psi_per_lsv_junction[.id == input_list[4]] > 0.1) %>%
  #dplyr::filter((mean_psi_per_lsv_junction[.id == input_list[4]] - mean_psi_per_lsv_junction[.id == input_list[3]] > 0.05)) %>% # |
                   #mean_psi_per_lsv_junction[.id == input_list[4]] - mean_psi_per_lsv_junction[.id == input_list[1]] > 0.05)) %>%
  mutate(color_gene_name = ifelse(mean_psi_per_lsv_junction[.id == input_list[4]] - mean_psi_per_lsv_junction[.id == input_list[3]] > 0.05, "rescued", "not_rescued")) %>%
  select(c(1,2,5,17,19)) %>%
  pivot_wider(names_from = ".id", values_from = "mean_psi_per_lsv_junction")
  #dplyr::filter(.id == "Control_Control")

input_list <- c("ctrl_ctrl",
                "ctrl_UPF1", 
                "TDP43_ctrl", 
                "TDP43_UPF1")
big_data_upf1_filtered <- big_data_upf1 %>%
  group_by(lsv_junc) %>%
  dplyr::filter(n() == 4) %>% #???
  dplyr::filter(de_novo_junctions == 1) %>%
  dplyr::filter(mean_psi_per_lsv_junction[.id == input_list[1]] < 0.05) %>%
  dplyr::filter(mean_psi_per_lsv_junction[.id == input_list[2]] < 0.05) %>%
  #dplyr::filter(mean_psi_per_lsv_junction[.id == input_list[3]] < 0.05) %>% 
  dplyr::filter(mean_psi_per_lsv_junction[.id == input_list[4]] > 0.1) %>%
  #dplyr::filter((mean_psi_per_lsv_junction[.id == input_list[4]] - mean_psi_per_lsv_junction[.id == input_list[3]] > 0.05)) %>% # |
                   #mean_psi_per_lsv_junction[.id == input_list[4]] - mean_psi_per_lsv_junction[.id == input_list[1]] > 0.05)) %>%
  mutate(color_gene_name = ifelse(mean_psi_per_lsv_junction[.id == input_list[4]] - mean_psi_per_lsv_junction[.id == input_list[3]] > 0.05, "rescued", "not_rescued")) %>%
  select(c(1,2,5,17,19)) %>%
  pivot_wider(names_from = ".id", values_from = "mean_psi_per_lsv_junction")
  #dplyr::filter(.id == "ctrl_ctrl")

big_data_combo <- inner_join(big_data_chx_filtered, big_data_upf1_filtered, by = c("lsv_junc","gene_name", "paste_into_igv_junction")) %>%
  filter(!grepl("\\.", gene_name)) %>%
  group_by(paste_into_igv_junction) %>% 
  slice_max(Control_TDP43KD)


library(viridis)
big_data_combo_matrix <- big_data_combo %>%
  filter(gene_name != "CORO7-PAM16") %>%
  #filter(Control_TDP43KD < 0.95) %>%
  select(c(3,5:8,10:13)) %>%
  column_to_rownames("paste_into_igv_junction") %>%
  as.matrix()
my_gene_col <- data.frame(rescue_chx = ifelse(big_data_combo_matrix[,4] - big_data_combo_matrix[,3] > 0.05, "rescued", "not_rescued"),
                          rescue_upf = ifelse(big_data_combo_matrix[,8] - big_data_combo_matrix[,7] > 0.05, "rescued", "not_rescued"))
table(my_gene_col)

my_colour = list(
  rescue_chx = c("rescued" = "#D53E4F", "not_rescued" = "#3288BD"),
  rescue_upf = c("rescued" = "#D53E4F", "not_rescued" = "#3288BD")
)

#mat_breaks <- seq(min(big_data_combo_matrix), max(big_data_combo_matrix), length.out = 10)

pheatmap::pheatmap(big_data_combo_matrix, 
                   color = colorRampPalette(rev(brewer.pal(n = 9, name = "Spectral")))(100),
                   annotation_colors = my_colour,
                   show_rownames = F,
                   #color = inferno(8),
                   #breaks = mat_breaks,
                   annotation_row = my_gene_col)



big_delta_upf1_filtered2 <- big_delta_upf1_filtered %>%
  #dplyr::filter(.id == "CtrlCtrl-TDP43Ctrl") %>%
  dplyr::filter(paste_into_igv_junction %in% big_delta_chx_filtered$paste_into_igv_junction)
big_delta_chx_filtered2 <- big_delta_chx_filtered %>%
  #dplyr::filter(.id == "ControlControl-ControlTDP43KD") %>%
  dplyr::filter(paste_into_igv_junction %in% big_delta_upf1_filtered$paste_into_igv_junction)

ggplot(data = big_delta_chx_filtered2, aes(x = reorder(paste_into_igv_junction, color_gene_name), y = 0, fill = color_gene_name)) +
  geom_tile() +
   geom_tile(data = big_delta_upf1_filtered2, aes(x = paste_into_igv_junction, y = -1, fill = color_gene_name)) +
  theme_classic() +
  colorspace::scale_fill_continuous_diverging() +
  theme(axis.text.x = element_text(angle = 90))


####
big_splicing_merged <- inner_join()
ploss <- big_splicing_merged %>%
  dplyr::filter(.id %in% c(input_list[3],input_list[4])) %>%
  ggplot(mapping = aes(x = .id, y = mean_psi_per_lsv_junction)) +
  facet_wrap(facets = vars(color_gene_name)) +
  geom_point(aes(color = color_gene_name, alpha = alpha_gene_name, group = lsv_junc), show.legend = F) + 
  geom_line(aes(color = color_gene_name, alpha = alpha_gene_name, group = lsv_junc), show.legend = F) +
  #geom_text_repel(aes(x = .id, y = mean_psi_per_lsv_junction, label = label_junction), point.padding = 0.3,
  #                nudge_y = 0.2, min.segment.length = 0.5, box.padding  = 2, max.overlaps = Inf, size=4, show.legend = F) +
  scale_color_manual(values = c(#"#ABDDA4", 
    "#3288BD",
    "#D53E4F")) +
  #scale_alpha_manual(values = c(0.02,1)) +
  xlab("") +
  ylab("PSI") +
  #scale_x_discrete(labels = c("TDP43-KD", "TDP43-KD + NMD inhibition")) +
  theme_classic()





####events that are significantly changed in chx+tdp/tdp vs ctrl
####need to run majiq between chx+tdp vs ctrl (also deseq2!!)

####events that are significantly changed in chx+tdp vs tdp

  

