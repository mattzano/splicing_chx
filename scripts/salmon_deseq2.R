# The tximport pipeline. 
# So this code has these following input that can be uptaken:
# (1) necessary: the father directory where all the salmon files folder locate + metadata that has sample_name column containing the name of each salmon folder.
# (2) necessary: specify the species. mouse/human will lead to the usage of a online ensembl data. users can also input "none" and use "--gtf" to input their own reference columns. But it is super slow, any suggestion...?
# (3) optional: counts from abundance. There are 4 options. read ?tximport for more information about difference between these 4.
# (4) optional: transcript or not. default is FALSE and will output the gene-level summarization from salmon files. but if users want the output table to be at transcript level (just a output of the original data), this command can achieve.
# (5) optional: save directory. the default is user's working directory.
# (6) optional: this function can also use the imported data to perform a default deseq2 test. Default is FALSE.
# (7) optional: for deseq2, users can use --setcontrol and --setcomparegroup to specify which group of samples are comparing and which one is the control that other groups are comparing with.


#make option that intake the directory etc
#suppressPackageStartupMessages(require(optparse))


#option_list = list(
#  make_option(c("-s","--salmondirectory"),action="store",type="character",default=NA,help="the directory where all the salmon quant files locate"),
#  make_option(c("-m","--metadata"),action="store",type="character",default=NA,help="the directory where the metadata file locates"),
#  make_option(c("-t","--tx2gene"),action="store",type="character",default=NA,help="Tx2gene file for the GTF the data was aligned to"),
#  make_option(c("-o","--outputdir"),action="store",default=getwd(), help="set the directory you want to save the results"),
#  make_option(c("-g","--column_name"),action="store", default="NA",help="Column name containing the metadata of baseline and contrasts that will be used"),
#  make_option(c("-b","--baseline"),action="store",default="baseline", help="What the control is called in the column name"),
#  make_option(c("-c","--contrast"),action="store",default="baseline", help="What the contrast is called in the column name"),
#  make_option(c("-n","--controls_name"),action="store",default="baseline", help="What the control should be called in the output files"),
#  make_option(c("-t","--contrast_name"),action="store",default="baseline", help="What the contrast should be called in the output files")
#)
#opt = parse_args(OptionParser(option_list=option_list))



#directory
#salmon_quant_directory = opt$salmondirectory
#outputdir = opt$outputdir
#metadata_dir = opt$metadata
#tx2gene = opt$tx2gene
#column_name = opt$column_name
#baseline = opt$baseline
#contrast = opt$contrast
#controls_name = opt$controls_name
#contrast_name = opt$contrast_name
#column_name = opt$column_name

outputdir = "~/Documents/Github/splicing_chx/results/deseq2/"
#metadata_dir = "~/Documents/Github/tdp43_concentration/data/metadata_dz_curves.csv" ###make this
tx2gene = "~/Desktop/rbp_bed/gencode.v40.tx2gene.csv"
column_name = "condition"
baseline = 0
contrast = 1
controls_name = "Control"
contrast_name = "TDP43KD"
#column_name = opt$column_name

output_path=paste0(outputdir,controls_name,"-",contrast_name,".")

# ============================ section 1: import data ===========================

tx2gene <- data.table::fread(tx2gene,header=FALSE)
tx2gene <- tx2gene[,c(1:2)]
colnames(tx2gene) = c("TXNAME", "GENEID")

#(1) First read in the metadata. if only a subset of the files are used, the opt$pattern option will be taken.

salmon_deseq2 <- function(mutation, sample_dir, metadata_dir) {
  
metadata2 = metadata_dir %>% 
  #dplyr::select(sample, !!as.symbol(column_name)) %>% 
  mutate(comparison_condition = case_when(!!as.symbol(column_name) == baseline ~ 'baseline',
                                          !!as.symbol(column_name) == contrast ~ 'contrast',
                                          TRUE ~ NA_character_)) %>% 
  dplyr::filter(!is.na(comparison_condition))

##

metadata <- metadata2 %>%
  #dplyr::mutate(group = gsub('[[:digit:]]+', '', sample)) #%>%
  dplyr::filter(treatment %in% mutation)
print(metadata)
    
    
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
                               design = ~ comparison_condition) 



# 'Note that the tximport-to-DESeq2 approach uses estimated gene counts from the transcript abundance quantifiers, but not normalized counts' -- <Deseq2 vignette> (just a note - Deseq() wraps the normalization step inside)
# perform the Deseq function
dds = DESeq(dds)

# Now, extract the result and named them by their contrast group
results_table <<- results(dds) %>%
  as.data.frame() %>%
  tibble::rownames_to_column('Geneid') %>%
  mutate(ensgene = gsub("\\..*", "", Geneid)) %>%
  left_join(annotables::grch38[c(1,3)])
  
# Now, extract the DESeq2 normed counts
normed_counts <<- counts(dds, normalized = TRUE) %>%
  as.data.frame() %>%
  tibble::rownames_to_column('Geneid') %>%
  mutate(ensgene = gsub("\\..*", "", Geneid)) %>%
  left_join(annotables::grch38[c(1,3)])
summary(normed_counts)
#return(results_table)

}

#write.csv(results_table,"DESeq2_results.csv")
#write.csv(normed_counts,"DESeq2_normalized_counts.csv")

