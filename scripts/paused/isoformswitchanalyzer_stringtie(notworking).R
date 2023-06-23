library(IsoformSwitchAnalyzeR)
salmon_quant_directory = "~/Desktop/stringtie"
metadata <- c("CTRL_ctrl_1","CTRL_ctrl_2","CTRL_ctrl_3","CTRL_ctrl_4",
              "DOX_ctrl_1", "DOX_ctrl_2", "DOX_ctrl_3", "DOX_ctrl_4")
files = unique(file.path(salmon_quant_directory,metadata,paste0("chx_tdp_",metadata,".assembled.gtf.tmap"))) 
names(files) = unique(metadata)

#import counts and abundance
stringTieQuant <- importIsoformExpression(sampleVector = files)
#create design table with conditions
myDesign <- data.frame(sampleID = colnames(stringTieQuant$abundance)[-1],
                       ###generalize this !!
                       condition = rep(c("CTRL_ctrl","DOX_ctrl"), each = 4))

#make switchAnalyzeRlist
aSwitchList <- importRdata(isoformCountMatrix   = stringTieQuant$counts,
                           isoformRepExpression = stringTieQuant$abundance,
                           designMatrix         = myDesign,
                           isoformExonAnnoation = "~/Desktop/gencode.v40.annotation.gtf",
                           isoformNtFasta       = "~/Desktop/gencode.v40.transcripts.fa")

#prefiltering
exampleSwitchListFiltered <- preFilter(switchAnalyzeRlist = aSwitchList)
#testing for isoform swithces
exampleSwitchListAnalyzed <- isoformSwitchTestDEXSeq(switchAnalyzeRlist = exampleSwitchListFiltered)
extractSwitchSummary(exampleSwitchListAnalyzed)

#extract sequences
"orfAnalysis" %in% names( exampleSwitchListFiltered )
exampleSwitchListAnalyzed <- extractSequence(exampleSwitchListAnalyzed,
                                             #removeLongAAseq = TRUE, alsoSplitFastaFile = TRUE,                                          
                                             pathToOutput = "~/Documents/GitHub/splicing_chx/results")

####external analysis bit - how do I include this to the pipeline?

#analyzeCPC2()
exampleSwitchListAnalyzed <- analyzeCPC2(switchAnalyzeRlist   = exampleSwitchListAnalyzed,
                                         pathToCPC2resultFile = "~/Documents/GitHub/splicing_chx/data/result_cpc2.txt",
                                         removeNoncodinORFs   = TRUE)   # because ORF was predicted de novo)

#analyzePFAM() - stopped for now (slow and needs excessive tweaking of data)
#exampleSwitchListAnalyzed <- analyzePFAM(switchAnalyzeRlist   = exampleSwitchListAnalyzed,
#  pathToPFAMresultFile = "~/Documents/GitHub/splicing_chx/data/pfam_results.txt",
#  showProgress=FALSE)

#analyzeSignalP()
exampleSwitchListAnalyzed <- analyzeSignalP(switchAnalyzeRlist       = exampleSwitchListAnalyzed,
                                            pathToSignalPresultFile  = "~/Documents/GitHub/splicing_chx/data/output_protein_type.txt")

#analyzeIUPred2A()
exampleSwitchListAnalyzed <- analyzeIUPred2A(switchAnalyzeRlist        = exampleSwitchListAnalyzed,
                                             pathToIUPred2AresultFile = "~/Documents/GitHub/splicing_chx/data/iupred2_result.txt",
                                             showProgress = FALSE)

#print brief report
exampleSwitchListAnalyzed

#predicting alternative splicing
exampleSwitchListAnalyzed <- analyzeAlternativeSplicing(switchAnalyzeRlist = exampleSwitchListAnalyzed, quiet=TRUE)
table( exampleSwitchListAnalyzed$AlternativeSplicingAnalysis$IR )

#change here
consequencesOfInterest=c(
  'intron_retention',
  'coding_potential',
  'ORF_seq_similarity',
  'NMD_status',
  #'domains_identified',
  'IDR_identified',
  'IDR_type',
  'signal_peptide_identified'
)

exampleSwitchListAnalyzed <- analyzeSwitchConsequences(exampleSwitchListAnalyzed,
                                                       consequencesToAnalyze = consequencesOfInterest, 
                                                       dIFcutoff = 0.1,
                                                       showProgress=FALSE
)

extractSwitchSummary(
  exampleSwitchListAnalyzed,
  filterForConsequences = TRUE
)
extractTopSwitches(exampleSwitchListAnalyzed, filterForConsequences = TRUE, n=100)

switchPlot(
  exampleSwitchListAnalyzed,
  gene='SNAP25',
  condition1 = 'CTRL_ctrl',
  condition2 = 'DOX_ctrl',
  localTheme = theme_bw(base_size = 13) # making text sightly larger for vignette
)

library(ggrepel)
filtered_data <- exampleSwitchListAnalyzed$isoformFeatures %>%
  dplyr::filter(abs(dIF) > 0.1 & isoform_switch_q_value < 0.00005 )
ggplot(data=exampleSwitchListAnalyzed$isoformFeatures, aes(x=dIF, y=-log10(isoform_switch_q_value))) +
  geom_point(
    aes( color=abs(dIF) > 0.1 & isoform_switch_q_value < 0.05 ), # default cutoff
    size=1
  ) +
  geom_text_repel(aes(label = gene_name), data = filtered_data) +
  geom_hline(yintercept = -log10(0.05), linetype='dashed') + # default cutoff
  geom_vline(xintercept = c(-0.1, 0.1), linetype='dashed') + # default cutoff
  facet_wrap( ~ condition_2) +
  #facet_grid(condition_1 ~ condition_2) + # alternative to facet_wrap if you have overlapping conditions
  scale_color_manual('Signficant\nIsoform Switch', values = c('black','red')) +
  labs(x='dIF', y='-Log10 ( Isoform Switch Q Value )') +
  theme_bw()

ggplot(data=exampleSwitchListAnalyzed$isoformFeatures, aes(x=gene_log2_fold_change, y=dIF)) +
  geom_point(
    aes( color=abs(dIF) > 0.1 & isoform_switch_q_value < 0.05 ), # default cutoff
    size=1
  ) + 
  facet_wrap(~ condition_2) +
  #facet_grid(condition_1 ~ condition_2) + # alternative to facet_wrap if you have overlapping conditions
  geom_hline(yintercept = 0, linetype='dashed') +
  geom_vline(xintercept = 0, linetype='dashed') +
  scale_color_manual('Signficant\nIsoform Switch', values = c('black','red')) +
  labs(x='Gene log2 fold change', y='dIF') +
  theme_bw()

#extractConsequenceSummary()
extractConsequenceSummary(
  exampleSwitchListAnalyzed,
  consequencesToAnalyze='all',
  plotGenes = FALSE,           # enables analysis of genes (instead of isoforms)
  asFractionTotal = FALSE      # enables analysis of fraction of significant features
)
#extractConsequenceEnrichment() and extractSplicingEnrichment()
extractConsequenceEnrichment(
  exampleSwitchListAnalyzed,
  consequencesToAnalyze='all',
  analysisOppositeConsequence = TRUE,
  returnResult = FALSE # if TRUE returns a data.frame with the summary statistics
)
#extractConsequenceEnrichmentComparison() - only when 3 or more conditions are analyzed
extractConsequenceEnrichmentComparison(
  exampleSwitchListAnalyzed,
  #consequencesToAnalyze=c('domains_identified','intron_retention','coding_potential'),
  analysisOppositeConsequence = TRUE,
  returnResult = FALSE # if TRUE returns a data.frame with the summary statistics
)
#extractConsequenceGenomeWide()
extractConsequenceGenomeWide(
  exampleSwitchListAnalyzed,
  featureToExtract = 'all',                 # all isoforms stored in the switchAnalyzeRlist
  #splicingToAnalyze = c('A3','MES','ATSS'), # Splice types significantly enriched in COAD
  plot=TRUE,
  returnResult=FALSE  # Preventing the summary statistics to be returned as a data.frame
)

#alternative splicing
extractSplicingSummary(
  exampleSwitchListAnalyzed,
  asFractionTotal = FALSE,
  plotGenes=FALSE
)

splicingEnrichment <- extractSplicingEnrichment(
  exampleSwitchListAnalyzed,
  splicingToAnalyze='all',
  returnResult=TRUE,
  returnSummary=TRUE
)

extractSplicingGenomeWide(
  exampleSwitchListAnalyzed,
  featureToExtract = 'all',                 # all isoforms stored in the switchAnalyzeRlist
  #splicingToAnalyze = c('A3','MES','ATSS'), # Splice types significantly enriched in COAD
  plot=TRUE,
  returnResult=FALSE  # Preventing the summary statistics to be returned as a data.frame
)