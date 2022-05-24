library(dplyr)
library(purrr)
library(tibble)
library(ggplot2)
library(ggVennDiagram)
library(ggvenn)
library(eulerr)

####euler diagram
bigdelta_filter <- big_data %>%
  group_by(lsv_junc) %>%
  dplyr::filter(.id == "Control_Control" | .id == "Control_TDP43KD") %>%
  dplyr::filter(n() == 2) %>%
  dplyr::filter(mean_psi_per_lsv_junction[.id == "Control_Control"] < 0.05) %>%
  dplyr::filter(mean_psi_per_lsv_junction[.id == "Control_TDP43KD"] > 0.1)

bigdelta_filter2 <- big_data %>%
  group_by(lsv_junc) %>%
  dplyr::filter(.id == "Control_Control" | .id == "Cycloheximide_TDP43KD") %>%
  dplyr::filter(n() == 2) %>%
  dplyr::filter(mean_psi_per_lsv_junction[.id == "Control_Control"] < 0.05) %>%
  dplyr::filter(mean_psi_per_lsv_junction[.id == "Cycloheximide_TDP43KD"] > 0.1)

bigdelta_filter22 <- bigdelta_filter2 %>%
  group_by(lsv_junc) %>%
  dplyr::filter(.id == "Cycloheximide_TDP43KD")

bigdelta_list <-list(TDP_only = bigdelta_filter$lsv_junc, #$lsv_junc[bigdelta_filter$.id == "Cycloheximide_TDP43KD"],
                     TDP_NMD = bigdelta_filter2$lsv_junc) #$lsv_junc[bigdelta_filter$.id == "Control_TDP43KD"])  

bigdelta_filter_merged <- rbind(bigdelta_filter,bigdelta_filter22)

#ggVennDiagram(bigdelta_list) +
#  ggplot2::scale_fill_continuous("blue", "red")
ggvenn(bigdelta_list, fill_color = c("#0073C2FF", "#EFC000FF"))

mat <- euler(c("TDP_only" = 1422, "TDP_NMD" = 7062,
               "TDP_only&TDP_NMD" = 1803))
plot(mat,
     quantities = list(type = c("percent", "counts")),
     fills = c("#0073C2FF", "#EFC000FF"))

#slope plot - why again????
big_data_filtereds <- big_data %>%
  group_by(lsv_junc) %>%
  dplyr::filter(n() == 4) %>%#
  dplyr::filter(mean_psi_per_lsv_junction[.id == "Control_Control"] < 0.05) %>%
  dplyr::filter(mean_psi_per_lsv_junction[.id == "Cycloheximide_Control"] < 0.05) %>%
  #dplyr::filter(mean_psi_per_lsv_junction[.id == "Control_TDP43KD"] > 0.05) %>%  ###remove this
  dplyr::filter(mean_psi_per_lsv_junction[.id == "Cycloheximide_TDP43KD"] > 0.05) %>%
  dplyr::filter((mean_psi_per_lsv_junction[.id == "Control_TDP43KD"] - mean_psi_per_lsv_junction[.id == "Control_Control"] > 0.05 |
                   mean_psi_per_lsv_junction[.id == "Cycloheximide_TDP43KD"] - mean_psi_per_lsv_junction[.id == "Control_Control"] > 0.05)) %>%
  mutate(color_gene_name = as.character(ifelse(mean_psi_per_lsv_junction[.id == "Cycloheximide_TDP43KD"] - mean_psi_per_lsv_junction[.id == "Control_TDP43KD"] > 0.05, 
                                               "Delta PSI > 0.05", "Delta PSI < 0.05"))) %>%
  mutate(alpha_gene_name = as.character(ifelse(gene_name %in% gene_list, 1, 0))) %>%#c("chr19:17641556-17642414","chr8:79611214-79616822"),"chr20:63439708-63444659","chr10:3099819-3101365","chr2:241668985-241670726","chr19:7169831-7170537"), 
  mutate(label_junction = case_when(.id =="Control_TDP43KD" & gene_name %in% gene_list ~ gene_name, T ~ ""))
table(big_data_filtereds$color_gene_name)

big_data_filtereds %>%
  dplyr::filter(.id %in% c("Control_TDP43KD","Cycloheximide_TDP43KD")) %>%
  ggplot(mapping = aes(x = .id, y = mean_psi_per_lsv_junction)) +
  facet_wrap(facets = vars(color_gene_name)) +
  geom_point(aes(color = color_gene_name, group = lsv_junc), alpha = 0.02, show.legend = F) + 
  geom_line(aes(color = color_gene_name, group = lsv_junc), alpha = 0.02, show.legend = F) +
  #geom_text_repel(aes(x = .id, y = mean_psi_per_lsv_junction, label = label_junction), point.padding = 0.3,
  #                nudge_y = 0.2, min.segment.length = 0.5, box.padding  = 2, max.overlaps = Inf, size=4, show.legend = F) +
  scale_color_manual(values = c(#"#ABDDA4", 
    "#3288BD",
    "#D53E4F")) +
  #scale_alpha_manual(values = c(0.02,1)) +
  xlab("") +
  ylab("PSI") +
  scale_x_discrete(labels = c("TDP43-KD", "TDP43-KD+CHX")) +
  scale_y_continuous(breaks = seq(0, 1, 0.1)) +
  theme_classic()

####use deseq2 to filter events - maybe use one from another dataset such as i3N or even proteomics (bottom of script)
my.dds <- run_standard_deseq("data/feature_count", 
                             base_grep = "CTRL",
                             contrast_grep = "DOX",  
                             grep_pattern = "ctrl",
                             baseName = "Untreated",
                             contrastName = 'TDP43KD')
results <- my.dds$results_table

my.dds2 <- run_standard_deseq("data/feature_count/chx_ctrl", 
                             base_grep = "_CTRL_ctrl_",
                             contrast_grep = "_DOX_chx_",  
                             grep_pattern = "",
                             baseName = "Untreated",
                             contrastName = 'TDP43KD')
results2 <- my.dds2$results_table

results_merged <- merge(results, results2, by  = "gene_name")

results_edit <- results_merged %>%
  #dplyr::filter(padj.x < 0.05) %>%
  group_by(gene_name) %>%
  mutate(delta_ctrl = log2FoldChange.x) %>%
  mutate(delta_chx = log2FoldChange.y) %>%
  mutate(color_gene_name = as.character(ifelse(delta_ctrl > 0 & delta_chx > 0, 1,
                                               ifelse(delta_ctrl > 0 & delta_chx < 0, 2,
                                                      ifelse(delta_ctrl < 0 & delta_chx > 0, 3, 4))))) %>%
  mutate(alpha_gene_name = as.character(ifelse(abs(delta_ctrl) < 2.5, # & abs(delta_chx) < 2.5, 
                                               1, 2))) %>%
  mutate(label_junction = case_when(alpha_gene_name == 2 ~ gene_name, T ~ ""))

results_edit %>%
  #dplyr::filter(gene_name == "CYFIP2" | gene_name == "UNC13A" | gene_name == "STMN2") %>%
  ggplot(aes(x = delta_ctrl, y = delta_chx)) +
  geom_text_repel(aes(label = label_junction), show.legend = F, max.overlaps = 100) +
  geom_point(aes(color = color_gene_name, alpha = alpha_gene_name), show.legend = F) +
  #geom_smooth(method='lm', formula= y~x) +
  scale_color_manual(values = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3")) +
  scale_alpha_manual(values = c(0.1, 0.8)) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  #eom_abline(slope = 1, intercept = 0) +
  coord_fixed() +
  theme_classic()



####similar plot but with psi instead of expression - gives lots of information but i don't like the dataviz
big_data_filterereds <- big_data %>%
  group_by(lsv_junc) %>%
  dplyr::filter(n() == 4) %>%
  dplyr::filter(mean_psi_per_lsv_junction[.id == "Control_Control"] < 0.05) %>%
  dplyr::filter(mean_psi_per_lsv_junction[.id == "Cycloheximide_Control"] < 0.05) %>%
  ###dplyr::filter(mean_psi_per_lsv_junction[.id == "Control_TDP43KD"] > 0.05) %>%
  dplyr::filter(mean_psi_per_lsv_junction[.id == "Cycloheximide_TDP43KD"] > 0.05) %>%
  dplyr::filter((mean_psi_per_lsv_junction[.id == "Control_TDP43KD"] - mean_psi_per_lsv_junction[.id == "Control_Control"] > 0.05 |
                   mean_psi_per_lsv_junction[.id == "Cycloheximide_TDP43KD"] - mean_psi_per_lsv_junction[.id == "Control_Control"] > 0.05)) %>%
  #mutate(color_gene_name = as.character(ifelse(mean_psi_per_lsv_junction[.id == "Cycloheximide_TDP43KD"] - mean_psi_per_lsv_junction[.id == "Control_TDP43KD"] > 0.05, 
  #                                             "Delta PSI > 0.05", "Delta PSI < 0.05")))
  mutate(delta_ctrl = mean_psi_per_lsv_junction[.id == "Control_TDP43KD"] - mean_psi_per_lsv_junction[.id == "Control_Control"]) %>%
  mutate(delta_chx = mean_psi_per_lsv_junction[.id == "Cycloheximide_TDP43KD"] - mean_psi_per_lsv_junction[.id == "Control_Control"]) %>%
  mutate(color_gene_name = as.character(ifelse(delta_ctrl > 0.1 & delta_chx > 0.1, 1,
                                               ifelse(delta_ctrl > 0.1 & delta_chx < 0.1, 2,
                                                      ifelse(delta_ctrl < 0.1 & delta_chx > 0.1, 3, 4))))) %>%
  mutate(alpha_gene_name = as.character(ifelse(delta_ctrl < 0.1 & delta_chx > 0.1, 1, 2))) %>%
  mutate(label_junction = case_when(alpha_gene_name == 1 ~ gene_name, T ~ ""))
table(big_data_filterereds$color_gene_name)


big_data_filterereds %>%
  dplyr::filter(.id %in% c("Control_Control")) %>%
  ggplot(aes(x = delta_ctrl, y = delta_chx)) +
  geom_text_repel(aes(label = label_junction), show.legend = F, max.overlaps = 100) +
  geom_point(aes(color = color_gene_name, alpha = alpha_gene_name), show.legend = F) +
  #geom_smooth(method='lm', formula= y~x) +
  scale_color_manual(values = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3")) +
  scale_alpha_manual(values = c(1, 0.1)) +
  #geom_vline(xintercept = 0) +
  #geom_hline(yintercept = 0) +
  #geom_abline(slope = 1, intercept = 0) +
  coord_fixed() +
  theme_classic()



#get data from patients or other DE data (proteomics?)
bbbb <- big_data_filterereds %>%
  dplyr::filter(.id %in% c("Control_Control")) %>%
  dplyr::filter(label_junction != "")

bbbb_to_bed <- bbbb[,c(10,13,14,2,11)]
write.table(bbbb_to_bed, "darkcryptome.bed", quote = F, sep = "\t", row.names = F, col.names = F)

dark_cryptome_counts <- read_csv("data/dark_cryptome_counts.csv") ####update?
dark_cryptome_counts %>%
  dplyr::filter(gene == "SYNE1") %>%
  ggplot(aes(x = disease, y = spliced_reads, color = tdp_path)) +
  geom_point(position = position_jitter(width = 0.35, height = 0.1), show.legend = F) +
  facet_wrap(facets = vars(region)) +
  theme_classic()
table(dark_cryptome_counts$disease[dark_cryptome_counts$region == "Spinal_Cord"])
####"C20orf194" "MTMR3" "SYNE1" "CYFIP2"
####"DLGAP3" "ERI3" "USP6NL" "DIP2C" "MEIS2" "AC069368.3" "PLEKHO2" "RPTOR" "CEP131" "UNC13B"
#####"MAST1" "CELSR3" "HTT" "FAM114A2" "CYFIP2"  "IGF2R" "KDM1B" "TAF1"

#tab <- read.table("data/nygc_dark_cryptomeaggregated.psi.tsv", sep = "\t")
#tab <- separate(tab, V1, into = c("PSI","sample_id","name","individual","region","tissue","tissue_clean","disease","disease_full","age",
#                                  "onset","mutations","pathology"), sep = ",")
#nygc_dark_cryptomeaggregated_psi_split <- separate(tab, name, into = c("gene", "category", "seq"), sep = "\\|")
#nygc_dark_cryptomeaggregated_psi_split %>%
  #dplyr::filter(disease != "NA" & disease != "Other") %>%
  #dplyr::filter(region != "NA" & region != "Other") %>%
#  dplyr::filter(gene == "CYFIP2") %>%
#  ggplot(aes(x = disease, y = PSI, color = category)) +
#  geom_point(position = position_dodge(width = 0.4), show.legend = T) +
#  facet_wrap(facets = vars(region)) +
#  theme_classic()
#######"ANKMY2" "DENND1B" "KDM1B" "ZNF583" "FAM114A2" "DCHS1"


###use proteomics data
library(readxl)
proteomics <- read.table("data/20200618_TDP43_KD_TotalProteome_DDA.csv", sep = ",", header = T)

#proteomics2 <- 
proteomics_small <- proteomics %>%
  dplyr::mutate(gene_color = as.character(ifelse(Log2FC > 0, 1, 2))) %>%
  dplyr::mutate(gene_alpha = as.character(ifelse((Log2FC > 1 | Log2FC < -1) & Minus_log10_Pvalue > 1, 1, 2))) %>%
  dplyr::mutate(gene_label = case_when(gene_alpha == 1  ~ Gene.names, T ~ "")) %>%
  dplyr::select("Gene.names", "Lead.gene.name", "Minus_log10_Pvalue", "Log2FC", "gene_color", "gene_alpha", "gene_label")
  
proteomics_small %>%
  #dplyr::filter(!is.nan(Log2FC)) %>%
  ggplot(aes(x = Log2FC, y = Minus_log10_Pvalue, color = gene_color, alpha = gene_alpha)) +
  geom_point() +
  geom_text_repel(aes(label = gene_label), color = "black", alpha = 1, show.legend = F, max.overlaps = 100) +
  scale_color_manual(values = c("#E41A1C", "#377EB8")) +
  scale_alpha_manual(values = c(1, 0.1)) +
  theme_classic()
 
darkcryptome_small <- big_data %>%
  group_by(lsv_junc) %>%
  dplyr::filter(n() == 4) %>%
  dplyr::filter(mean_psi_per_lsv_junction[.id == "Control_Control"] < 0.05 & 
                  mean_psi_per_lsv_junction[.id == "Cycloheximide_Control"] < 0.05 & 
                  mean_psi_per_lsv_junction[.id == "Cycloheximide_TDP43KD"] > 0.05) %>%
  dplyr::filter((mean_psi_per_lsv_junction[.id == "Control_TDP43KD"] - mean_psi_per_lsv_junction[.id == "Control_Control"] > 0.05 |
                   mean_psi_per_lsv_junction[.id == "Cycloheximide_TDP43KD"] - mean_psi_per_lsv_junction[.id == "Control_Control"] > 0.05)) %>%
  mutate(delta_ctrl = mean_psi_per_lsv_junction[.id == "Control_TDP43KD"] - mean_psi_per_lsv_junction[.id == "Control_Control"]) %>%
  mutate(delta_chx = mean_psi_per_lsv_junction[.id == "Cycloheximide_TDP43KD"] - mean_psi_per_lsv_junction[.id == "Control_Control"]) %>%
  mutate(color_gene_name = as.character(ifelse(delta_ctrl > 0.1 & delta_chx > 0.1, 1,
                                               ifelse(delta_ctrl > 0.1 & delta_chx < 0.1, 2,
                                                      ifelse(delta_ctrl < 0.1 & delta_chx > 0.1, 3, 4))))) %>%
  mutate(alpha_gene_name = as.character(ifelse(color_gene_name == 3, 1, 2))) %>%
  mutate(label_junction = case_when(alpha_gene_name == 1 ~ gene_name, T ~ "")) %>%
  dplyr::filter(.id == "Control_Control") %>%
  dplyr::select(".id", "gene_name", "mean_psi_per_lsv_junction", "de_novo_junctions", "lsv_junc", "paste_into_igv_junction", "exon_type", "delta_ctrl", "delta_chx", "color_gene_name",
                "alpha_gene_name", "label_junction")

prot_crypt <- left_join(proteomics_small, darkcryptome_small, by = c("Lead.gene.name" = "gene_name"))

prot_crypt_plot <- prot_crypt %>%
  dplyr::mutate(label = case_when(alpha_gene_name == 1 & gene_alpha == 1 ~ Lead.gene.name, T ~ ""))

prot_crypt_plot %>%
  ggplot(aes(x = Log2FC, y = delta_chx-delta_ctrl, color = color_gene_name, alpha = gene_alpha)) +
  geom_point() + 
  geom_text_repel(aes(label = label), color = "black", alpha = 1, show.legend = F, max.overlaps = 1000) +
  #geom_hline(yintercept = 0) +
  scale_color_manual(values = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3"),
                     labels = c("DeltaPSI > 0.1 in both", "DeltaPSI < 0.1 only in TDP43KD", "Delta PSI > 0.1 only with TDP43KD+CHX", "DeltaPSI < 0.1 in both")) +
  scale_alpha_manual(values = c(1, 0.1),
                     labels = c("Significantly up/downregualted in proteomics", "Not significant")) +
  xlab("Log2FC in proteomics") + ylab("DeltaPSI TDP43KD+CHX - DeltaPSI TDP43KD") +
  theme_classic()


