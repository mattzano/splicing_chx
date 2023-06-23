postar_bed <- fread("~/Desktop/rbp_bed/human_RBP_binding_sites_sorted.bed", 
                    col.names = c("seqnames", "start", "end", "dataset", "score", "strand", "QC"))

sy5y_tdp_bed_1 <- fread("~/Desktop/rbp_bed/tardbp-shsy5y-1-20210701-mh_mapped_to_genome_single_peaks.bed", 
                      col.names = c("seqnames", "start", "end", "dataset", "score", "strand"))
sy5y_tdp_bed_1 <- sy5y_tdp_bed_1 %>%
  mutate(QC = "-") %>%
  mutate(dataset = "zSy5yTDP")

sy5y_tdp_bed_2 <- fread("~/Desktop/rbp_bed/tardbp-shsy5y-2-20210701-mh_mapped_to_genome_single_peaks.bed", 
                        col.names = c("seqnames", "start", "end", "dataset", "score", "strand"))
sy5y_tdp_bed_2 <- sy5y_tdp_bed_2 %>%
  mutate(QC = "-") %>%
  mutate(dataset = "zSy5yTDP")

postar_bed <- rbind(postar_bed, sy5y_tdp_bed_1, sy5y_tdp_bed_2)


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



gene_target <- c("CYFIP2")

cds_parent = as_tibble(cds_regions[cds_regions$transcript_id %in% princ$ensembl_transcript_id[princ$external_gene_name == gene_target]])
exons_parent = as_tibble(exons_regions[exons_regions$transcript_id %in% princ$ensembl_transcript_id[princ$external_gene_name == gene_target]])

parent_cryptic_delta  <- big_delta %>%  #change to input_splicing when adding to the Rmd
  dplyr::filter(gene_name == gene_target & .id == "CycloheximideControl-CycloheximideTDP43KD" & mean_dpsi_per_lsv_junction > 0 & probability_changing > 0.9) %>%
  group_by(paste_into_igv_junction) %>%
  dplyr::filter(mean_dpsi_per_lsv_junction > 0.1 & base_mean_psi < 0.05)

parent_postar_bed <- postar_bed %>%
  dplyr::filter(seqnames %in% exons_parent$seqnames & start > min(exons_parent$start) & end < max(exons_parent$end)) %>%
  #dplyr::filter(seqnames %in% parent_cryptic_delta$seqnames & start > min(parent_cryptic_delta$start)-1000 & end < max(parent_cryptic_delta$end)+1000) %>%
  mutate(RBP = paste0(".",word(dataset, 1, sep = "_"))) %>%
  mutate(colour_gene = as.character(ifelse(QC == "-", 1, 0))) %>%
  group_by(RBP) %>% dplyr::filter(n() >= 20 | RBP == ".zSy5yTDP") 

#if (nrow(parent_cryptic_delta) > 0) {# && nrow(exons_parent) > 0  && nrow(cds_parent) > 0) {
  
tabbb <- to_intron(cds_parent)
plotz <- ggplot(aes(xstart = start, xend = end, y = gene_target), data = exons_parent) +
  geom_range(data = exons_parent, fill = "white", height = 0.2) +
  geom_range(data = cds_parent, fill = "black", height = 0.4) +
  geom_intron(data = to_intron(cds_parent), aes(strand = strand)) +
  geom_junction(aes(color = junc_cat), data = parent_cryptic_delta, junction.orientation = "top", show.legend = T) +
  geom_range(aes(y=RBP, fill = colour_gene), data = parent_postar_bed, height = 0.3, show.legend = F) +
    #geom_range(aes(y=0.5), data = parent_sy5y_tdp_bed, color = "#A37EB8", fill = "#A37EB8", height = 0.3) +
  #ggforce::facet_zoom(xlim = c(min(parent_cryptic_delta$start)-500, max(parent_cryptic_delta$end)+500)) +
  ylab("") +
  scale_y_discrete(expand = c(0,2)) +
  scale_x_continuous() +
  theme_classic2() +
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank())
  
plot(plotz)