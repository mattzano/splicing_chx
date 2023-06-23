library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(tidyverse)
library(data.table)
library(GenomicFeatures)
library(plyranges)
## 'cds_by_tx' must be a GRangesList object obtained with cdsBy(txdb, by="tx")
addCdsPhase <- function(cds_by_tx)
{
    cds_phase <- pc(rep(IntegerList(0), length(cds_by_tx)),
                    heads((3L - (cumsum(width(cds_by_tx)) %% 3L)) %% 3L, n=-1L))
    unlisted_cds_by_tx <- unlist(cds_by_tx, use.names=FALSE)
    mcols(unlisted_cds_by_tx)$cds_phase <- unlist(cds_phase, use.names=FALSE)
    relist(unlisted_cds_by_tx, cds_by_tx)
}

format_grange_for_gtf = function(grange,novel_event,gene_transcript){
    
    if(is.character(novel_event)){
        string_cryp = gsub(":|-","_",novel_event,1,-3)
    }else{
        string_cryp = gsub(":|-","_",stringr::str_sub(novel_event$novel_exon,1,-3))
        
    }
    
    new_transcript_id = glue::glue("{gene_transcript}_{string_cryp}")
    new_transcript_id = glue::glue("{gene_transcript}_{string_cryp}")
    
    grange$transcript_id = new_transcript_id
    
    grange = GeneStructureTools::reorderExonNumbers(grange)
    
    grange$exon_rank = grange$exon_number
    grange$original_transcript = gene_transcript
    ##the gene
    
    gene_id_this_gene = human_genes[transcript == gene_transcript,unique(gene_id.y)]
    grange$tx_biotype = human_genes[transcript == gene_transcript,unique(TXBIOTYPE)]
    grange$gene_name = human_genes[transcript == gene_transcript,unique(gene_name)]
    
    grange$gene_id = gene_id_this_gene
    
    names(grange) = ifelse(!is.na(grange$exon_name),grange$exon_name,grange$novel_exon)
    
    return(grange)
}

my_orf = function (transcripts, BSgenome = NULL, returnLongestOnly = TRUE, 
                   allFrames = FALSE,  
                   uORFs = FALSE) 
{
    if (allFrames == TRUE) {
        returnLongestOnly = FALSE
        longest = 1
    }
    
    transcripts$exon_number <- as.numeric(transcripts$exon_number)
    order <- order(transcripts$transcript_id, transcripts$exon_number)
    transcripts <- transcripts[order]
    transcripts$seq <- as.character(Biostrings::getSeq(BSgenome, 
                                                       transcripts))
    seqCat <- aggregate(seq ~ transcript_id, mcols(transcripts), 
                        function(x) (paste(x, collapse = "")))
    ids <- as.character(seqCat$transcript_id)
    
    seqCat <- seqCat$seq
    rm <- which(grepl("N", seqCat))
    if (length(rm) > 0) {
        seqCat <- seqCat[-rm]
        removeId <- ids[rm]
        ids <- ids[-rm]
        transcripts <- transcripts[-which(transcripts$transcript_id %in% 
                                              removeId)]
    }
    seqCat <- c(seqCat, stringr::str_sub(seqCat, 2), stringr::str_sub(seqCat, 
                                                                      3))
    frames <- rep(c(1, 2, 3), each = length(ids))
    ids <- c(ids, ids, ids)
    orf <- suppressWarnings(unlist(lapply(seqCat, function(x) as.character(Biostrings::translate(Biostrings::DNAString(x))))))
    orfDF <- data.frame(id = ids, aa_sequence = orf, frame = frames, 
                        stringsAsFactors = FALSE)
    orfDF$seq_length <- nchar(orfDF$aa_sequence)
    orfDF$seq_length_nt <- nchar(seqCat) + orfDF$frame - 1
    
    
    return(orfDF)
}




txdb = transcripts(TxDb.Hsapiens.UCSC.hg38.knownGene)
cryptic_regions = fread("/Users/annaleigh/Documents/GitHub/biomarker_with_sahba/data/concated_ce_with_transcripts.txt")

cryptic_regions_protein_coding = cryptic_regions %>% 
    mutate(V2 = ifelse(V6 == "-",V2 + 1,V2)) %>% 
    mutate(V2 = ifelse(V6 == "+",V2 + 1,V2)) %>% 
    left_join(annotables::grch38_tx2gene, 
              by = c("V8" = "enstxp")) %>% 
    left_join(annotables::grch38 %>% select(ensgene,biotype)) %>% 
    left_join(annotables::grch38 %>% select(symbol,biotype),
              by = c("V7" = "symbol")) %>% 
    mutate(biotype_full = case_when(is.na(biotype.y) ~ biotype.x,
                                    is.na(biotype.x) ~ biotype.y,
                                    is.na(biotype.x) & is.na(biotype.y) ~ NA_character_,
                               TRUE ~ biotype.y)) %>% 
    select(-biotype.x,-biotype.y) %>% 
    filter(biotype_full == 'protein_coding')

cryptic_regions_protein_coding = cryptic_regions_protein_coding %>% 
    makeGRangesFromDataFrame(keep.extra.columns = TRUE,
                             seqnames.field = 'V1',
                            start.field = 'V2',
                            end.field = 'V3',
                            strand.field = 'V6')

cryptic_regions_protein_coding$transcript_id = cryptic_regions_protein_coding$V8

cds_regions = cdsBy(TxDb.Hsapiens.UCSC.hg38.knownGene, "tx",use.names = TRUE)
cds_regions = unlist(cds_regions)

cds_regions$transcript_id = gsub("\\..*", "", names(cds_regions))


rm(i)
peptides = c()
with_peptide = c()

for(i in 1:length(cryptic_regions_protein_coding)){

    gene_name = cryptic_regions_protein_coding[i]$V7
    print(gene_name)
    parent_transcript = cryptic_regions_protein_coding[i]$transcript_id
    cds_parent = cds_regions %>% filter(transcript_id == parent_transcript)
    
    if(length(cds_parent) == 0){
        next()
    }
    
    exon_one = cryptic_regions_protein_coding[i]
    
    new_model = sort(c(exon_one,cds_parent))
    
    new_model = GeneStructureTools::reorderExonNumbers(new_model)
    
    names(new_model) = NULL
    new_model = unlist(addCdsPhase(GRangesList(new_model)))
    
    cryptic_number = new_model %>% plyranges::filter(is.na(cds_id)) %>% 
        as.data.frame() %>% pull(exon_number)
    
    cryptic_up_down = new_model %>% filter(exon_number %in% c(cryptic_number - 1, 
                                                              cryptic_number,cryptic_number + 1))
    
    cds_parent = GeneStructureTools::reorderExonNumbers(cds_parent)
    normal = my_orf(cds_parent,BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38) %>% 
        dplyr::slice(1) %>% pull(aa_sequence)

    # tmp = GeneStructureTools::getOrfs(new_model,BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38,
    #                                   returnLongestOnly = FALSE, 
    #                                   allFrames = TRUE,
    #                                   uORFs = TRUE)
    tmp = my_orf(new_model,BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38)
    tmp$gene_id = gene_name
    # tmp$nmd_prob <- notNMD::predictNMD(tmp, "prob")
    ce_aa = tmp %>% dplyr::slice(1) %>% pull(aa_sequence)
    if(gene_name == "AARS"){
        break()
    }
    
    if(str_count(ce_aa,"\\*") == 1){
        aln = msa::msa(Biostrings::AAStringSet(list(Biostrings::AAString(normal), Biostrings::AAString(ce_aa))))
        ce_pos = which(strsplit(msa::msaConsensusSequence(aln)[[1]],"")[[1]] == "?")
        ce_peptide = paste0(strsplit(ce_aa,"")[[1]][ce_pos],collapse = "")
        print(gene_name)
        print(ce_peptide)
        peptides = c(peptides,ce_peptide)
        with_peptide = c(with_peptide,gene_name)
    }
    
}




exon_seq <- Biostrings::translate(Biostrings::getSeq(BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38, 
                                                     cryptic_regions_protein_coding[1]))
df <- data.frame(seqnames=seqnames(new_model),
                 starts=start(cryptic_regions_protein_coding)-1,
                 ends=end(new_model),
                 names=c(rep(".", length(new_model))),
                 scores=c(rep(".", length(new_model))),
                 strands=strand(new_model))

write.table(df, file="foo.bed", quote=F, sep="\t", row.names=F, col.names=F)



