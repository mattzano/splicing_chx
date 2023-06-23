#1- Download bigwig file of conservation track (via biomaRt/AnnotationHub in R). #Import to R via rtracklayer::import.bw()
#2- Overlap your regions of interest with the the conservation object #(e.g. GenomicRanges join_overlap_* functions)

#phast <- phastCons100way.UCSC.hg38
#class(phast)

###necessary files
big_delta
postar_bed <- fread("~/Desktop/rbp_bed/human_RBP_binding_sites_sorted.bed", 
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




gene_list = "UNC13A"
conservation(gene_list)

conservation <- function(gene_list) {
gene_list_ensgene <- annotables::grch38 %>%
  dplyr::filter(symbol == gene_list) %>%
  pull(ensgene)

cds_parent = as_tibble(cds_regions[cds_regions$transcript_id %in% princ$ensembl_transcript_id[princ$external_gene_name == gene_list]])

exons_parent = as_tibble(exons_regions[exons_regions$transcript_id %in% princ$ensembl_transcript_id[princ$external_gene_name == gene_list]])

parent_cryptic_delta  <- big_delta %>%  #change to input_splicing when adding to the Rmd
  dplyr::filter(gene_name == gene_list & 
                (.id == "ControlControl-ControlTDP43KD" | .id == "CycloheximideControl-CycloheximideTDP43KD") & 
                mean_dpsi_per_lsv_junction > 0 & probability_changing > 0.9) %>%
  group_by(paste_into_igv_junction) %>%
  dplyr::filter(mean_dpsi_per_lsv_junction > 0.1 & base_mean_psi < 0.05)

parent_postar_bed <- postar_bed %>%
  dplyr::filter(seqnames %in% exons_parent$seqnames & start > min(exons_parent$start) & end < max(exons_parent$end)) %>%
  dplyr::filter(seqnames %in% parent_cryptic_delta$seqnames & start > min(parent_cryptic_delta$start)-500 & end < max(parent_cryptic_delta$end)+500) %>%
  mutate(RBP = paste0(".",word(dataset, 1, sep = "_")))

poldip3 <- subset(CNCscore_GRanges, overlap_any_gene_v92_name %in% gene_list_ensgene)
poldip3_df <- data.frame(poldip3@ranges@start, poldip3$CDTS, poldip3$mean_phastCons20way, poldip3$Log2Ratio_CDTS_cons)
names(poldip3_df) <- c("location", "CDTS", "phyloP", "ratio")

p <- ggplot() +
  geom_area(aes(x = location, y = CDTS), 
            data = poldip3_df, fill = "red", alpha = 0.2) +
  geom_area(aes(x = location, y = phyloP), 
            data = poldip3_df, fill = "blue", alpha = 0.5) +
  geom_area(aes(x = location, y = ratio), 
            data = poldip3_df, alpha = 0.5) +
  #geom_area(aes(x = location, y = value), 
  #          data = cons, alpha = 0.5, fill = "red") +
  #geom_line(aes(x = start, y = CDTS), 
  #          data = cdts_poldip3, color = "red", alpha = 0.5) +
  geom_range(aes(xstart = start, xend = end, y = 0), 
             data = exons_parent, fill = "white", height = 0.2) +
  geom_range(aes(xstart = start, xend = end, y = 0), 
             data = cds_parent, fill = "black", height = 0.4) +
  geom_intron(aes(xstart = start, xend = end, y = 0, strand = unique(exons_parent$strand)), 
              data = to_intron(cds_parent)) +
  ylab("phyloP (blue), CDTS (red), and CNC (grey) scores") +
  xlab(paste0("Genomic location of ", gene_list, " gene in ", unique(exons_parent$seqnames))) +
  theme_classic() +
  theme(axis.line.y = element_blank(), 
        axis.ticks.y = element_blank())

if (nrow(parent_cryptic_delta) > 0) {
  p +
    geom_junction(aes(xstart = start, xend = end, y = 0), 
                  data = parent_cryptic_delta, show.legend = F, junction.orientation = "top") +
    geom_range(aes(xstart = start, xend = end, y=0, fill = RBP, colour = RBP), 
               data = parent_postar_bed, height = 0.1, show.legend = T) +
    ggforce::facet_zoom(xlim = c(min(parent_cryptic_delta$start)-500, max(parent_cryptic_delta$end)+500))
} else {
  p
}
#score_tab <- gscores(phast, GRanges(seqnames = unique(exons_parent$seqnames),
#                                    ranges= IRanges(start=min(exons_parent$start):max(exons_parent$end),width = 1)))
#score_df <- data.frame(score_tab@ranges@start,score_tab$default)
#names(score_df) <- c("location", "value")

#cons <- read.table("~/Downloads/unc13a.txt", header = F,skip = 9)
#names(cons) <- c("location", "value")
#rolled10 <- rollmean(cons$value,10)
#cons$rolled <- c(rep(0,5), rolled10, rep(0,4))
#cons_binned <- cons %>%
#  group_by(G=trunc(2:(n()+1)/10)) %>% 
#  summarise(value=mean(value),location=min(location)) %>% 
#  dplyr::select(-G)
#CDTS9 <- read.table("~/Downloads/chr9.txt", header = T) ####changehere
#cdts_poldip3 <- CDTS9 %>%
  #mutate(start = as.numeric(start)) %>%
  #mutate(end = as.numeric(end)) %>%
#  dplyr::filter(start >= 84668375) %>%
#  dplyr::filter(end <= 84753850)

#ggplot() +
#  geom_area(aes(x = location, y = value), 
#            data = score_df, alpha = 0.5) +
  #geom_area(aes(x = location, y = value), 
  #          data = cons, alpha = 0.5, fill = "red") +
  #geom_line(aes(x = start, y = CDTS), 
  #          data = cdts_poldip3, color = "red", alpha = 0.5) +
#  geom_range(aes(xstart = start, xend = end, y = 0), 
#             data = exons_parent, fill = "white", height = 0.2) +
#  geom_range(aes(xstart = start, xend = end, y = 0), 
#             data = cds_parent, fill = "black", height = 0.4) +
#  geom_intron(aes(xstart = start, xend = end, y = 0, strand = unique(exons_parent$strand)), 
#              data = to_intron(cds_parent)) +
  #geom_junction(aes(xstart = start, xend = end, y = 0), 
  #             data = parent_cryptic_delta, show.legend = F, junction.orientation = "top") +
#  geom_range(aes(xstart = start, xend = end, y=0, fill = RBP, colour = RBP), 
#             data = parent_postar_bed, height = 0.1, show.legend = T) +
  #ggforce::facet_zoom(xlim = c(min(parent_cryptic_delta$start)-500, max(parent_cryptic_delta$end)+500)) +
#  ylab("phyloP score") +
#  theme_classic() +
#  theme(axis.line.y = element_blank(), 
#        axis.ticks.y = element_blank())
}
