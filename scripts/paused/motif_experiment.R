###motif tdp43 binding 

#start_gr <- GRanges(seqnames = "chr1",
#                    ranges= IRanges(start=10000:11000,width = 1001))
#end_gr <- GRanges(seqnames = "chr1",
#                  ranges= IRanges(start=110000:111000,width = 1001))
#chr1_vect <- seq(from = 1, to = 248956422, by = 100000)

#sequences <- GRanges(seqnames = "chr1",
#                    ranges= IRanges(start=chr1_vect,width = 10000))
#df <- data.frame(sequences@seqnames,sequences@ranges)
#names(df)[1] <- "seqnames"
#names(df)[6] <- "strand"
#df$sequence_chr1 <- as.character(getSeq(g, sequences))
#sequence_chr1 <- data.frame(sequence = getSeq(g, sequences))
#motif_detected <- df %>% dplyr::filter(grepl("TG.ATG",sequence_chr1))


#cnc_genes <- #CNCscore_GRanges[ CNCscore_GRanges$overlap_any_gene_v92_name != NA ]
#  subset(CNCscore_GRanges, !is.na(overlap_any_gene_v92_name))
#cnc_chr1 <- subset(CNCscore_GRanges, seqnames == "chr1")
#cnc_chr1_genes <- subset(cnc_chr1, !is.na(overlap_any_gene_v92_name))
#df <- data.frame(cnc_chr1_genes@seqnames,cnc_chr1_genes@ranges)
#names(df)[1] <- "seqnames"


width <- 10
motif = c("TG") #,"GT")
window = 3
o <- org.Hs.eg.db

g <- BSgenome.Hsapiens.UCSC.hg38
tx <- TxDb.Hsapiens.UCSC.hg38.knownGene
tx_chr1 <- genes(tx)
tx_chr2 <- subset(tx_chr1, seqnames %in% c("chr1"))#,"chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10",                  
                                           #"chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20",                  
                                           #"chr21","chr22","chrX","chrY"))
df <- data.frame(tx_chr2@seqnames,tx_chr2@ranges,tx_chr2@strand,tx_chr2$gene_id)
names(df)[1] <- "seqnames"
names(df)[6] <- "strand"
names(df)[7] <- "gene_id"
df$sequence_chr1 <- as.character(getSeq(g, tx_chr2))
motif_detected <- df %>% dplyr::filter(grepl("TG.ATG",sequence_chr1))

###find "GU" and "UG" motifs around the splicing sites "GURAGU" and "...."
motif_detected$seq_chr1 <- str_split(motif_detected$sequence_chr1, "TG.ATG")
motif_detected$seq_chr2 <- str_count(motif_detected$sequence_chr1, "TG")
motif_detected$seq_chr3 <- str_count(motif_detected$sequence_chr1, "GT")
motif_detected$length <- motif_detected$end-motif_detected$start
motif_detected$ratio <- (motif_detected$seq_chr2 + motif_detected$seq_chr3) / motif_detected$length


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


#select_seqs <- function(sequences, motif = "TG"){
#motif_detected <- dplyr::filter(sequences, grep(motif, sequence))
#return(motif_detected)
#}


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
  
  #start_result <- select_seqs(start_seqs, motif = motif)
  motif_detected <- start_seqs %>% dplyr::filter(grepl(motif,sequence))
  
  
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


####
#get list of genes - with metadata!
#found genomic location
#filter for splicing sites
#find tdp43 binding + / - 1000 bp
#spliceAI




####splice ai around those motifs

###see whether enrichment for neuronal genes