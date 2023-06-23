library(tidyverse)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(patchwork)

find_motif <- function(sequences, width, motif = "TG"){
  matrix <- matrix(0, nrow = length(sequences), ncol = str_length(sequences[1]))
  
  motif_detected <- str_locate_all(sequences, motif)
  for(i in 1:length(motif_detected)){
    for(j in motif_detected[i]){
      matrix[i,j] <- 1
    }
  }
  return(colSums(matrix)/length(sequences))
}

make_gr <- function(seqnames, positions, width){
  
  yo <- GRanges(seqnames = seqnames,
                ranges = IRanges(start = positions-width,
                                 end = positions+width))
  
  return(yo)
}

make_motif_df <- function(g, df, width = 1000, motif = "TG..TG", window = 3){
  # This function takes a df with the columns "chr", junc_start" and "junc_end" CHANGED THAT
  
  start_gr <- make_gr(df$seqnames, df$start, width)
  end_gr <- make_gr(df$seqnames, df$end, width)
  
  start_seqs <- data.frame(sequence = getSeq(g, start_gr))
  end_seqs <- data.frame(sequence = getSeq(g, end_gr))
  
  
  start_result <- find_motif(start_seqs$sequence, width = width, motif = motif)
  plot(start_result)
  
  end_result <- find_motif(end_seqs$sequence, width = width, motif = motif)
  plot(end_result)
  
  
  df2 <- data.frame(donor = start_result, acceptor = end_result) %>%
    mutate(pos = -width:width) %>%
    pivot_longer(cols = c("donor", "acceptor")) %>%
    group_by(name) %>%
    mutate(smooth = zoo::rollmean(value, k= window, na.pad=T))
  
  plot <- ggplot(df2, aes(pos, smooth, colour = name)) +
    geom_line() +
    facet_wrap(~name)
  print(plot)
  
  return(df2)
}

width <- 1000
motif = c("TG") #,"GT")
window = 3
g <- BSgenome.Hsapiens.UCSC.hg38


######need to get cryptic exons not junctions
cryptic_df<- big_delta %>%  #change to input_splicing when adding to the Rmd
  dplyr::filter(.id == "ControlControl-ControlTDP43KD") %>%
  dplyr::filter(probability_changing > 0.9)
  #group_by(paste_into_igv_junction) %>%
  #dplyr::filter(mean_dpsi_per_lsv_junction > 0.1 & base_mean_psi < 0.05) 
chr19:17642414-17642541

positive <- cryptic_df %>%
  dplyr::filter(mean_dpsi_per_lsv_junction > 0.1) #%>%
  #dplyr::filter(base_mean_psi < 0.05) %>%
  #dplyr::filter(probability_changing > 0.9) #%>%
  #dplyr::filter(strand == "-")

positive_df <- make_motif_df(g, positive, window = 3)
df = data.frame(c("seqnames","chr19"),
                c("start", 17642414), 
                c("end",17642541))
df <- row_to_names(df,1)
df[2] <- as.numeric(df[2])

negative <- cryptic_df %>%
  dplyr::filter(mean_dpsi_per_lsv_junction < -0.1) #%>%
  #dplyr::filter(probability_changing > 0.9) #%>%
  #dplyr::filter(strand == "-")

negative_df <- make_motif_df(g, negative, window = 3)


### background 
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
annotated <- data.frame(intronsByTranscript(txdb, use.names=TRUE))

sub <- annotated %>% 
  dplyr::filter(str_length(seqnames) <= 4) %>% # remove weird things
  sample_frac(0.1) %>%
  #dplyr::filter(strand == "-") %>% ########################## ?
  dplyr::select(seqnames, start, end)

background_df <- make_motif_df(g, sub, window = 6)


pos_compare <- positive_df %>%
  mutate(type = "positive_delta") %>%
  bind_rows(background_df %>% mutate(type = "background"))

pos_compare_plot <- ggplot(pos_compare, aes(x = pos, y = smooth, colour = type)) +
  geom_line() +
  facet_grid(~factor(name, levels = c("donor", "acceptor"))) +
  ggtitle("Prevalence of UGNNUG around splice sites with increased PSI upon TDP43 KD") +
  theme_bw()

neg_compare <- negative_df %>%
  mutate(type = "negative_delta") %>%
  bind_rows(background_df %>% mutate(type = "background"))

neg_compare_plot <- ggplot(neg_compare, aes(x = pos, y = smooth, colour = type)) +
  geom_line() +
  facet_grid(~factor(name, levels = c("donor", "acceptor"))) +
  ggtitle("Prevalence of UGNNUG around splice sites with decreased PSI upon TDP43 KD") +
  theme_bw()

combined_plot <- pos_compare_plot/neg_compare_plot
combined_plot
#ggsave(plot = combined_plot, "combined_UGNNUG_very_wide.pdf")
#ggsave(plot = combined_plot, "combined_UGNNUG_very_wide.png")


###spliceai -I input.vcf -O output.vcf -R genome.fa -A grch37
