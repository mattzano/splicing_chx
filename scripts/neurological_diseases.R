###panel gene list

input_disease <- list.files(path = here::here("~/Documents/GitHub/splicing_chx/data/panels"))
datalist_disease <- list()

for (i in input_disease) {
  #input_tsv <- paste(i, "tsv", sep = ".")
  input_disease_tsv <- read.table(file.path("~/Documents/GitHub/splicing_chx/data/panels", i), sep = "\t", header = T, fill = TRUE)
  input_disease_tsv$Model_Of_Inheritance <- gsub("\\,.*", "", input_disease_tsv$Model_Of_Inheritance)
  input_disease_clean <- input_disease_tsv[,c(3,5,6,7,8,22,23)]
  input_disease_clean$gene_name <- input_disease_clean$Gene.Symbol
  datalist_disease[[i]] <- input_disease_clean
}

big_disease <- rbindlist(datalist_disease, idcol = TRUE)
big_disease$.id <- factor(big_disease$.id, levels = input_disease)
big_disease <- big_disease %>%
  mutate(category = ifelse(Model_Of_Inheritance == "BIALLELIC", "BIALLELIC",
                         ifelse(Model_Of_Inheritance == "MONOALLELIC", "MONOALLELIC",
                                ifelse(Model_Of_Inheritance == "BOTH monoallelic and biallelic (but BIALLELIC mutations cause a more SEVERE disease form)" |
                                         Model_Of_Inheritance == "BOTH monoallelic and biallelic", "BOTH",
                                       ifelse(Model_Of_Inheritance == "X-LINKED: hemizygous mutation in males","X-LINKED","Other")))))


##for exomes
names(big_disease)[2] <- "symbol"  
match <- big_disease %>%
  left_join(annotables::grch38) %>%
  dplyr::filter(chr %in% c(1:22) | chr == "X") %>%
  dplyr::select(c(13:15))
#write.table(match, "neuro_disease.bed", row.names = F, col.names = F, quote = F)

###for chx
big_delta_clean <- big_delta_chx %>%
  dplyr::filter(.id == "ControlControl-ControlTDP43KD") %>%
  dplyr::filter(base_mean_psi < 0.05 & mean_dpsi_per_lsv_junction > 0.1 & probability_changing > 0.9)

big_delta_disease <- big_delta_clean %>%
  dplyr::filter(gene_name %in% big_disease$gene_name) 

big_delta_disease_labelled <- left_join(big_delta_disease, big_disease, by = "gene_name")
#big_delta_disease_$Model_Of_Inheritance <- gsub("\\,.*", "", big_delta_disease$Model_Of_Inheritance)

big_delta_disease_labelled %>%
  dplyr::filter(category == "BIALLELIC") %>%
  ggplot(aes(x = Level4, fill = Model_Of_Inheritance)) +
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90))


big_delta_disease_labelled2 <- big_delta_disease_labelled %>%
  mutate(gene_label = case_when((category == "BIALLELIC" | category == "X-LINKED" | category == "BOTH") ~ gene_name, T ~ ""))

big_delta_disease_labelled3 <- big_delta_disease_labelled2 %>% 
  distinct(Level3, gene_label, .keep_all = T)


big_delta_disease_labelled2 %>%
  ggplot(aes(x = Level3, y = mean_dpsi_per_lsv_junction)) +
  #facet_wrap(facets = vars(Level3), scales = "free_y") +
  geom_point(aes(color = category), position = position_dodge(width = 0.5)) +
  geom_text_repel(aes(label = gene_label), data = big_delta_disease_labelled3,  max.overlaps = Inf) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90))



new_cryptic_junctions <- read.table("~/Desktop/rbp_bed/new_cryptics.junctions.bed")

#new_cryptic_junctions$gene_name <- gsub("\\|.*","",new_cryptic_junctions$V4)

big_cryptic_disease <- new_cryptic_junctions %>%
  mutate(paste_into_igv_junction = paste0(V1,":",V2,"-",V3)) %>%
  separate(V4, into = c("gene_name", "junc_type", "count"), sep = "\\|") %>%
  dplyr::filter(gene_name %in% big_disease$gene_name) 

big_cryptic_disease_labelled <- left_join(big_cryptic_disease, big_disease, by = "gene_name")

big_cryptic_disease_labelled %>%
  dplyr::filter(Model_Of_Inheritance == "BIALLELIC" & junc_type != "annotated" & junc_type != "none" & junc_type != "ambig_gene") %>%
  ggplot(aes(x = Level4, fill = junc_type)) +
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90))

big_cryptic_disease_labelled2 <- big_cryptic_disease_labelled %>%
  mutate(gene_label = case_when(category == "BIALLELIC" | category == "X-LINKED" ~ gene_name, T ~ ""))

#big_cryptic_disease_labelled3 <- big_cryptic_disease_labelled2 %>% 
#  distinct(Level3, gene_label, .keep_all = T)

#label_data <- big_cryptic_disease_labelled2
# calculate the ANGLE of the labels
#number_of_bar <- nrow(label_data)
#angle <-  90 - 360 * (5-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
# calculate the alignment of labels: right or left
#label_data$hjust<-ifelse( angle < -90, 1, 0) # If I am on the left part of the plot, my labels have currently an angle < -90
#label_data$angle<-ifelse(angle < -90, angle+180, angle) # flip angle BY to make them readable


big_cryptic_disease_labelled2 %>%
  #dplyr::filter(Level3 == "Sleep disorders") %>%
  ggplot(aes(x = as.factor(Level3))) +
  #facet_wrap(facets = vars(Level3), scales = "free_y") +
  geom_bar(aes(fill = category), color = "black", position = position_dodge(width = 0.5)) +
  #geom_text_repel(aes(label = gene_label), data = big_cryptic_disease_labelled3,  max.overlaps = Inf) +
  theme_minimal() + theme(axis.text.x = element_text(angle = 90))
    #axis.text = element_blank(),
    #axis.title = element_blank(),
    #panel.grid = element_blank()) #,
    #plot.margin = unit(rep(-2,4), "cm")) #+
    #theme(axis.text.x = element_text(angle = 90)) +
#  ylim(-25,300) +
  #scale_y_continuous(trans = "sqrt", limits = c(-1000,300)) +
 # coord_polar(start = 0) +
#  geom_text(data=label_data, aes(x=as.factor(Level3), label = ..count.., hjust=hjust), stat = "count", 
#            color="black", fontface="bold",alpha=0.6, size=2.5, angle= label_data$angle, inherit.aes = FALSE )




changing_sign <- read.table("~/Desktop/rbp_bed/changing_sig.junctions.bed")
#changing_sign$gene_name <- gsub("\\|.*","",changing_sign$V4)
big_changing_disease <- changing_sign %>%
  mutate(paste_into_igv_junction = paste0(V1,":",V2,"-",V3)) %>%
  separate(V4, into = c("gene_name", "junc_type", "psi_ctrl", "psi_tdp"), sep = "\\|") %>%
  dplyr::mutate(psi_ctrl = as.numeric(psi_ctrl)) %>%
  dplyr::mutate(psi_tdp = as.numeric(psi_tdp)) %>%
  dplyr::mutate(delta_psi = psi_tdp-psi_ctrl) %>%
  dplyr::filter(gene_name %in% big_disease$gene_name)

big_changing_disease_labelled <- left_join(big_changing_disease, big_disease, by = "gene_name")

big_changing_disease_labelled %>%
  #dplyr::filter(Model_Of_Inheritance == "BIALLELIC" & junc_type != "annotated" & junc_type != "none" & junc_type != "ambig_gene") %>%
  ggplot(aes(x = Level4, fill = junc_type)) +
  facet_wrap(facets = vars(category), scales = "free_y") +
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90))

big_changing_disease_labelled2 <- big_changing_disease_labelled %>%
  mutate(gene_label = case_when(psi_ctrl < 0.05 & (category == "BIALLELIC" | category == "X-LINKED") ~ gene_name, T ~ ""))
  
big_changing_disease_labelled3 <- big_changing_disease_labelled2 %>% 
  distinct(Level3, gene_label, .keep_all = T)

big_changing_disease_labelled2 %>%
  #dplyr::filter(Level3 == "Sleep disorders") %>%
  ggplot(aes(x = Level3, y = delta_psi)) +
  #facet_wrap(facets = vars(Level3), scales = "free_y") +
  geom_point(aes(color = category), position = position_dodge(width = 0.5)) +
  geom_text_repel(aes(label = gene_label), data = big_changing_disease_labelled3,  max.overlaps = Inf) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90)) #+
  #coord_polar()

#table1 <- as.data.frame(table(big_disease$Level4))
#table2 <- as.data.frame(table(big_changing_disease_labelled$Level4))
#table3 <- left_join(table1,table2, by = "Var1")



