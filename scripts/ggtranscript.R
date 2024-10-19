#Biostrings::getSeq(BSgenome, transcripts)...
#Biostrings::translate...

#postar_bed <- fread("~/Desktop/rbp_bed/human_RBP_binding_sites_sorted.bed", 
#                    col.names = c("seqnames", "start", "end", "dataset", "score", "strand", "QC"))

#ensembl <- biomaRt::useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
#humandb <- biomaRt::getBM(attributes = c("external_gene_name", 
#                                         "ensembl_gene_id", "ensembl_transcript_id", "transcript_appris", 
#                                         "uniprotswissprot", "uniprotsptrembl", "uniprot_isoform"), 
#                          mart = ensembl)
#princ <- humandb[which(humandb$transcript_appris == "principal1" |
#                       humandb$transcript_appris == "principal2" | 
#                       humandb$transcript_appris == "principal3" | 
#                       humandb$transcript_appris == "principal4" |
#                       humandb$transcript_appris == "principal5"), ]

#cds_regions = cdsBy(TxDb.Hsapiens.UCSC.hg38.knownGene, "tx",use.names = TRUE)
#cds_regions = unlist(cds_regions)
#cds_regions$transcript_id = gsub("\\..*", "", names(cds_regions))

#exons_regions = exonsBy(TxDb.Hsapiens.UCSC.hg38.knownGene, "tx",use.names = TRUE)
#exons_regions = unlist(exons_regions)
#exons_regions$transcript_id = gsub("\\..*", "", names(exons_regions))

my_clean_reader = function(file){
  as.data.table(janitor::clean_names(fread(file)))
}
my_name_fixer = function(tbl){
  colnames(tbl) = gsub(colnames(tbl)[[9]],"psi",colnames(tbl))
  return(tbl)
}


transcript_bind_plot <- function(gene_target) {
 #gene_target = "HDGFL2"
  cds_parent = as_tibble(cds_regions[cds_regions$transcript_id %in% princ$ensembl_transcript_id[princ$external_gene_name == gene_target]])
  exons_parent = as_tibble(exons_regions[exons_regions$transcript_id %in% princ$ensembl_transcript_id[princ$external_gene_name == gene_target]])

# parent_cryptic <- psi_list_full %>%
#    dplyr::filter(seqnames %in% exons_parent$seqnames & start > min(exons_parent$start) & end < max(exons_parent$end)) %>%
#    #group_by(seqnames, start, end, strand_junction) %>%
#    dplyr::mutate(sample_category = paste0(word(SampleID, 2, sep = "_"), "-", word(SampleID, 3, sep = "_"))) %>%
#    dplyr::filter(sample_category == "dox-ctrl") %>% #| sample_category == "ctrl-ctrl") %>%
#    dplyr::filter(type != "annotated") %>%
#    dplyr::filter(psi > 0.1) %>%
#   dplyr::mutate(igv_junc = paste0(seqnames, ":", start, "-", end))
 
# parent_cryptic_nt <- psi_list_full %>%
#   dplyr::filter(seqnames %in% exons_parent$seqnames & start > min(exons_parent$start) & end < max(exons_parent$end)) %>%
#   #group_by(seqnames, start, end, strand_junction) %>%
#   dplyr::mutate(sample_category = paste0(word(SampleID, 2, sep = "_"), "-", word(SampleID, 3, sep = "_"))) %>%
#   dplyr::filter(sample_category == "ctrl-ctrl") %>%
#   dplyr::filter(type != "annotated") %>%
#   dplyr::filter(psi < 0.05) %>%
#   dplyr::mutate(igv_junc = paste0(seqnames, ":", start, "-", end))
 
# parent <- rbind(parent_cryptic, parent_cryptic_nt)

#cparent_cryptic_wider <- parent %>%
#   pivot_wider(id_cols = c(igv_junc, type),
#               names_from = SampleID,
#               values_from = psi)

#parent_cryptic_wider_filter <- parent_cryptic_wider %>%
#  mutate(mean_psi_dox = rowMeans(parent_cryptic_wider[2:5], na.rm = T)) %>%
#  mutate(mean_psi_ctrl = rowMeans(parent_cryptic_wider[6:9], na.rm = T)) %>%
#  mutate(delta_psi = mean_psi_dox - mean_psi_ctrl) %>%
#  dplyr::filter(!is.na((delta_psi)))

 parent_cryptic_delta  <- big_delta_chx %>%  #change to input_splicing when adding to the Rmd
   dplyr::filter(gene_name == gene_target & .id == "ControlControl-ControlTDP43KD" & 
                   mean_dpsi_per_lsv_junction > 0 & probability_changing > 0.9) %>%
   arrange(desc(mean_dpsi_per_lsv_junction)) %>% 
   distinct(gene_name, .keep_all = T)
   #group_by(paste_into_igv_junction) %>%
   #dplyr::filter(mean_dpsi_per_lsv_junction > 0.1 & base_mean_psi < 0.05)
 #parent_cryptic_delta <- parent_cryptic_delta[c(1,3),]   
 #datapp <- as.data.frame(c(parent_cryptic_delta[1,3], parent_cryptic_delta[2,4]))
  parent_postar_bed <- postar_bed %>%
    dplyr::filter(seqnames %in% exons_parent$seqnames & start > min(exons_parent$start) & end < max(exons_parent$end)) %>%
    dplyr::filter(seqnames %in% parent_cryptic_delta$seqnames & start > min(parent_cryptic_delta$start)-100 & end < max(parent_cryptic_delta$end)+100) %>%
    mutate(RBP = paste0(".",word(dataset, 1, sep = "_"))) %>%
    #mutate(colour_gene = as.character(ifelse(QC == "-", 1, 0))) %>%
    #group_by(RBP) %>% dplyr::filter(n() > 3) %>%
    mutate(colour_gene = ifelse(RBP == ".TARDBP", "red",
                                       "blue"))
    #dplyr::filter(RBP != ".TARDBP")
  
  #parent_sy5y_tdp_bed <- sy5y_tdp_bed %>%
  #  dplyr::filter(seqnames %in% exons_parent$seqnames & start > min(exons_parent$start) & end < max(exons_parent$end))

#datapp <- as.data.frame(c(parent_cryptic_delta[2,3], parent_cryptic_delta[1,4]))
  plotz <- ggplot(aes(xstart = start, xend = end, y = gene_target), data = exons_parent) +
      geom_range(data = exons_parent, fill = "white", height = 0.2) +
      geom_range(data = cds_parent, fill = "black", height = 0.4) +
      geom_intron(data = to_intron(cds_parent), aes(strand = strand)) +
      geom_junction(data = parent_cryptic_delta, show.legend = F, junction.orientation = "top") +
      geom_range(aes(y=RBP, fill = colour_gene, colour = colour_gene), data = parent_postar_bed, height = 0.3, show.legend = F) +
      ggforce::facet_zoom(xlim = c(min(parent_cryptic_delta$start)-100, max(parent_cryptic_delta$end)+100)) +
      #geom_range(data = datapp) +
       #geom_segment(aes(x = 79630000, xend = 79630000, y = 4.75, yend = 4.95), 
      #             arrow = arrow(type="closed", length = unit(5,"points"))) +
      ylab("") +
      scale_y_discrete(expand = c(0,2)) +
      scale_fill_brewer(palette = "Set2") +
      scale_colour_brewer(palette = "Set2") +
      theme_classic2() +
      theme(axis.line.y = element_blank(), axis.ticks.y = element_blank())
    
  if (nrow(parent_cryptic_delta) > 0 & nrow(parent_postar_bed) > 0) {
      plot(plotz)
  } else {
      print(paste0("MAJIQ found no cryptics in ", gene_target))
    
 }
}


#Biostrings::getSeq(BSgenome, transcripts)...
#Biostrings::translate...




#unc_cryptic  <- big_data %>%
#  dplyr::filter(gene_name == "UNC13A" & (.id == "Control_Control" | .id == "Control_TDP43KD")) %>%  #generalise gene name
#  group_by(paste_into_igv_junction, exon_type) %>%
#  dplyr::filter(mean_psi_per_lsv_junction[.id == "Control_TDP43KD"] > 0.1 &
#                  mean_psi_per_lsv_junction[.id == "Control_Control"] < 0.05) %>%
#  mutate(delta_psi = mean_psi_per_lsv_junction[.id == "Control_TDP43KD"] - mean_psi_per_lsv_junction[.id == "Control_Control"]) %>%
#  dplyr::filter(delta_psi > 0.05 & .id == "Control_TDP43KD")
#unc_conc <- fread("data/unc_conc.csv")
#names(unc_conc)[c(5,7,8)] <- c("seqnames", "start", "end")
#unc_merge <- AppendMe(c("unc_conc", "unc_cryptic_df"))