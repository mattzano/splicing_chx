library(magrittr)
library(dplyr)
library(ggplot2)
library(ggtranscript)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(GenomicRanges)
library(rtracklayer)
library(ensembldb)
library(stringr)
library(biomaRt)

#AppendMe <- function(dfNames) {
#  do.call(rbind, lapply(dfNames, function(x) {
#    cbind(get(x), source = x)
#  }))
#}
#Biostrings::getSeq(BSgenome, transcripts)...
#Biostrings::translate...

postar_bed <- fread("data/metagene/human_RBP_binding_sites_sorted.bed", 
                    col.names = c("seqnames", "start", "end", "dataset", "score", "strand", "QC"))

ensembl <- biomaRt::useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
humandb <- biomaRt::getBM(attributes = c("external_gene_name", 
                                         "ensembl_gene_id", "ensembl_transcript_id", "transcript_appris", 
                                         "uniprotswissprot", "uniprotsptrembl", "uniprot_isoform"), 
                          mart = ensembl)
princ <- humandb[which(humandb$transcript_appris == "principal1" |
                       humandb$transcript_appris == "principal2" | 
                       humandb$transcript_appris == "principal3" | 
                       humandb$transcript_appris == "principal4" |
                       humandb$transcript_appris == "principal5"), ]

cds_regions = cdsBy(TxDb.Hsapiens.UCSC.hg38.knownGene, "tx",use.names = TRUE)
cds_regions = unlist(cds_regions)
cds_regions$transcript_id = gsub("\\..*", "", names(cds_regions))

exons_regions = exonsBy(TxDb.Hsapiens.UCSC.hg38.knownGene, "tx",use.names = TRUE)
exons_regions = unlist(exons_regions)
exons_regions$transcript_id = gsub("\\..*", "", names(exons_regions))


transcript_bind_plot <- function(gene_target) {

  cds_parent = as_tibble(cds_regions[cds_regions$transcript_id %in% princ$ensembl_transcript_id[princ$external_gene_name == gene_target]])
  exons_parent = as_tibble(exons_regions[exons_regions$transcript_id %in% princ$ensembl_transcript_id[princ$external_gene_name == gene_target]])

  parent_cryptic  <- big_delta %>%  #change to input_splicing when adding to the Rmd
    dplyr::filter(gene_name == gene_target & .id == "ControlControl-ControlTDP43KD") %>%
    group_by(paste_into_igv_junction) %>%
    dplyr::filter(mean_dpsi_per_lsv_junction > 0.1 & base_mean_psi < 0.05)

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

  parent_postar_bed <- postar_bed %>%
    dplyr::filter(seqnames %in% parent_cryptic$seqnames & start > min(parent_cryptic$start)-1000 & end < max(parent_cryptic$end)+1000) %>%
    mutate(RBP = paste0(".",word(dataset, 1, sep = "_"))) %>%
    group_by(RBP) %>% dplyr::filter(n() > 3)

  if (nrow(parent_cryptic) > 0 & 
      nrow(parent_postar_bed) > 0 & 
      nrow(exons_parent) > 0  & 
      nrow(cds_parent) > 0) {
    exons_parent %>%
      ggplot(aes(xstart = start, xend = end, y = gene_target)) +
      geom_range(data = exons_parent, fill = "white", height = 0.2) +
      geom_range(data = cds_parent, fill = "black", height = 0.4) +
      geom_intron(data = to_intron(cds_parent), aes(strand = strand)) +
      geom_junction(data = parent_cryptic, color = "#E41A1C", show.legend = F, junction.y.max = 0.6) +
      geom_range(aes(y=RBP), data = parent_postar_bed, color = "#377EB8", fill = "#377EB8", height = 0.3) +
      ggforce::facet_zoom(xlim = c(min(parent_cryptic$start)-1000, max(parent_cryptic$end)+1000)) +
      ylab("") +
      scale_y_discrete() + #expand = c(0,1)) +
      theme_classic2() +
      theme(axis.line.y = element_blank(), axis.ticks.y = element_blank())
    }
}

#big_delta_filter <- big_delta %>%
#  dplyr::filter(.id == "ControlControl-ControlTDP43KD") %>%
#  group_by(paste_into_igv_junction) %>%
#  dplyr::filter(mean_dpsi_per_lsv_junction > 0.1 & base_mean_psi < 0.05)

#gene_list <- unique(big_delta_filter$gene_name[!big_delta_filter$gene_name %in% c("EVI5","IFT122","WASH6P","RP1-138B7.8","RP11-206L10.17",
#                                                                                  "RP11-101E3.5","RP11-286N22.8","RP4-583P15.15",
#                                                                                  "CH17-264B6.6","RP11-554A11.11","HERC2P3","SEPTIN7P2")]) 

#for (i in gene_list[88:91]) {
#  plot(transcript_bind_plot(i))
#}
