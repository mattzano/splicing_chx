my_clean_reader = function(file){
  as.data.table(janitor::clean_names(fread(file)))
}
my_name_fixer = function(tbl){
  colnames(tbl) = gsub(colnames(tbl)[[9]],"psi",colnames(tbl))
  return(tbl)
}
file_path = "data/sj_tabs/normalized_annotated/"
suffix = "_normalized_annotated.csv"
psi_files = list.files(file_path,
                       pattern = suffix,
                       full.names = TRUE)
psi_list_full = purrr::map(psi_files,my_clean_reader)
samp_ids = base::tolower(purrr::simplify(purrr::map(psi_files, basename)))
samp_ids = gsub(suffix,"",samp_ids)
##add on the name as an additional column
psi_list_full = purrr::map2(psi_list_full, samp_ids, ~cbind(.x, SampleID = .y))
psi_list_full = purrr::map(psi_list_full,my_name_fixer)

psi_list_full = data.table::rbindlist(psi_list_full)




psi_list_unc <- psi_list_full %>%
  dplyr::filter(seqnames %in% parent_cryptic$seqnames & start > min(exons_parent$start) & end < max(exons_parent$end)) %>%
  #group_by(seqnames, start, end, strand_junction) %>%
  dplyr::mutate(sample_category = paste0(word(SampleID, 2, sep = "_"), "-", word(SampleID, 3, sep = "_"))) %>%
  dplyr::filter(sample_category == "dox-ctrl" | sample_category == "ctrl-ctrl") %>%
  dplyr::filter(type != "annotated") # & type != "novel_exon_skip") #%>%
  #group_by(seqnames, start, end) %>%
  #dplyr::filter(n() >= 2) #%>%
  #dplyr::filter(sample_category == "ctrl-ctrl"

psi_list_unc_wider <- psi_list_unc %>%
  pivot_wider(id_cols = c(seqnames, start, end),
              names_from = SampleID,
              values_from = psi) #%>%

psi_list_unc_ctrl <- psi_list_unc %>%
  dplyr::filter(sample_category == "ctrl-ctrl" & psi < 0.05) %>%
  group_by(seqnames, start, end) %>%
  dplyr::filter(n() >= 2)

psi_list_unc_dox <-psi_list_unc %>%
  dplyr::filter(sample_category == "dox-ctrl" & psi > 0.1) %>%
  group_by(seqnames, start, end) %>%
  dplyr::filter(n() >= 2)
  #group_by(seqnames, start, end, sample_category) %>%
  #dplyr::mutate(mean_psi = mean(psi)) %>%
  #distinct(seqnames, start, end, sample_category, .keep_all = T)#%>%
  # %>%
  #unique()

#psi_list_unc <- psi_list_unc %>%
#  distinct(seqnames, start, end, sample_category, .keep_all = T)

parent_postar_bed_tdp <- parent_postar_bed %>%
  dplyr::filter(RBP == ".TARDBP")

psi_list_unc %>%
  #dplyr::filter(type != "annotated" & type != "novel_exon_skip") %>%
  ggplot(aes(xstart = start, xend = end, y = "UNC13A", color = type)) +
  geom_range(data = exons_parent, fill = "white", color = "black", height = 0.2) +
  geom_range(data = cds_parent, fill = "black", color = "black", height = 0.4) +
  #geom_junction() +
  geom_junction(data = psi_list_unc_ctrl, junction.orientation = "top") +
  geom_junction(data = psi_list_unc_dox, junction.orientation = "bottom") +
  #geom_intron(data = to_intron(cds_parent), aes(strand = strand)) +
  #geom_junction(data = parent_cryptic, color = "#E41A1C", show.legend = F, junction.y.max = 0.6) +
  geom_range(aes(y=RBP), data = parent_postar_bed_tdp, color = "#377EB8", fill = "#377EB8", height = 0.3) +
  #ggforce::facet_zoom(xlim = c(17630000, 17633000)) +
  ylab("") +
  scale_y_discrete(expand = c(0,1.5)) +
  theme_classic2() +
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank())
#plot(psi_list_unc)
