
carm <- read_excel("~/Documents/phd/research_lines/c9orf_mouse/analysis_c9mouse/data/Spinalcord_GR-PR-WT_Reanalysis_CM.xlsx")
#carm <- read_excel("~/Documents/phd/research_lines/c9orf_mouse/Cortex_GR-PR-WT_Reanalysis_Four PR replicates.xlsx")
carmm <- carm %>% dplyr::filter(Gene_names != "Nnt")



####volcano plot
carmm %>%
  mutate(log2fc = log2(FC_PR_WT)) %>%
  mutate(color_name = (ifelse(p_value_PR_WT > 0.05, "Not significant",
                             ifelse(log2fc > 0, "Upregulated", "Downregulated")))) %>%
  mutate(color_name = factor(color_name, levels = c("Not significant","Upregulated", "Downregulated"))) %>%
  mutate(junction_name = case_when(((log2fc > 1 | log2fc < -1) & p_value_PR_WT < 0.05) | Gene_names == "C9orf72" ~ Gene_names,
                                   T ~ "")) %>%
  mutate(graph_alpha = ifelse((log2fc > 1 | log2fc < -1) & p_value_PR_WT < 0.05, 1, 0.2)) %>%
  ggplot(aes(x =log2fc, y = `Log_p-value_PR_WT`)) +
  geom_point(aes(color = color_name), show.legend = T) +
  geom_text_repel(aes(label = junction_name),
                  point.padding = 0.3,
                  nudge_y = 0.2,
                  min.segment.length = 0,
                  box.padding  = 2,
                  size = 2,
                  #max.overlaps = 1000,
                  show.legend = F,
                  max.overlaps = Inf) +
  geom_hline(yintercept = -log10(0.05), size =0.1) +
  geom_vline(xintercept = 0, linetype="dotted") +
  scale_color_manual(values = c("#899DA4", "#3B9AB2", "#FAD510")) +
  #guides(alpha = FALSE, size = FALSE) +
  theme_classic() +
  labs(x = expression(paste("Protein abundance (lo", g[2], "FC)")),
       y = expression(paste("-lo", g[10], " (", italic("P"), " value)")),
       color = "",
       title = "Volcano plot from (PR)400 mouse proteomics") +
  theme(legend.position = "top",
        legend.justification='left') +
  #xlim(-2,2) +
  theme(text = element_text(size = 8)) +
  theme(legend.text=element_text(size=6))
#ggsave("~/Desktop/volcanoPR.pdf", width = 3.93, height = 4.245)

carmm %>%
  mutate(log2fc = log2(FC_GR_WT)) %>%
  mutate(color_name = (ifelse(p_value_GR_WT > 0.05, "Not significant",
                              ifelse(log2fc > 0, "Upregulated", "Downregulated")))) %>%
  mutate(color_name = factor(color_name, levels = c("Not significant","Upregulated", "Downregulated"))) %>%
  mutate(junction_name = case_when(((log2fc > 1 | log2fc < -1) & p_value_GR_WT < 0.05) | Gene_names == "C9orf72" ~ Gene_names,
                                   T ~ "")) %>%
  mutate(graph_alpha = ifelse((log2fc > 1 | log2fc < -1) & p_value_GR_WT < 0.05, 1, 0.2)) %>%
  ggplot(aes(x =log2fc, y = `Log_p-value_GR_WT`)) +
  geom_point(aes(color = color_name), show.legend = T) +
  geom_text_repel(aes(label = junction_name),
                  point.padding = 0.3,
                  nudge_y = 0.2,
                  min.segment.length = 0,
                  box.padding  = 2,
                  size = 2,
                  max.overlaps = Inf,
                  show.legend = F) +
  geom_hline(yintercept = -log10(0.05), size =0.1) +
  geom_vline(xintercept = 0, linetype="dotted") +
  scale_color_manual(values = c("#899DA4", "#3B9AB2", "#FAD510")) +
  #guides(alpha = FALSE, size = FALSE) +
  theme_classic() +
  labs(x = expression(paste("Protein abundance (lo", g[2], "FC)")),
       y = expression(paste("-lo", g[10], " (", italic("P"), " value)")),
       color = "",
       title = "Volcano plot from (GR)400 mouse proteomics") +
  theme(legend.position = "top",
        legend.justification='left') +
  #xlim(-2,2) +
  theme(text = element_text(size = 8)) +
  theme(legend.text=element_text(size=6))
#ggsave("~/Desktop/volcanoGR.pdf", width = 3.93, height = 4.245)



#create input for GO
cm_PR_up <- carmm %>% dplyr::filter(FC_PR_WT > 1 & p_value_PR_WT < 0.05) 
cm_PR_up <- cm_PR_up[order(-cm_PR_up$FC_PR_WT),]
cm_PR_down <- carmm %>% dplyr::filter(FC_PR_WT < 1 & p_value_PR_WT < 0.05)
cm_PR_down <- cm_PR_down[order(cm_PR_down$FC_PR_WT),]

cm_GR_up <- carmm %>% dplyr::filter(FC_GR_WT > 1 & p_value_GR_WT < 0.05)
cm_GR_up <- cm_GR_up[order(-cm_GR_up$FC_GR_WT),]
cm_GR_down <- carmm %>% dplyr::filter(FC_GR_WT < 1 & p_value_GR_WT < 0.05)
cm_GR_down <- cm_GR_down[order(cm_GR_down$FC_GR_WT),]

cm_all_up <- carmm %>% dplyr::filter(FC_PR_WT > 1 & p_value_PR_WT < 0.05 &
                                       FC_GR_WT > 1 & p_value_GR_WT < 0.05) 
cm_all_up <- cm_all_up[order(-cm_all_up$FC_GR_WT),]
cm_all_down <- carmm %>% dplyr::filter(FC_PR_WT < 1 & p_value_PR_WT < 0.05 &
                                       FC_GR_WT < 1 & p_value_GR_WT < 0.05) 
cm_all_down <- cm_all_down[order(cm_all_down$FC_GR_WT),]

#GO analysis
go_cm_all <- gprofiler2::gost(query = list(Upregulated = cm_all_up$Gene_names, 
                                           Downregulated = cm_all_down$Gene_names), 
                                sources = c("GO:BP", "GO:CC", "GO:MF"), organism = 'mmusculus', 
                                ordered_query = T, evcodes = TRUE)
go_cm_PR <- gprofiler2::gost(query = list(Upregulated = cm_PR_up$Gene_names,
                                          Downregulated = cm_PR_down$Gene_names), 
                             sources = c("GO:BP", "GO:CC", "GO:MF"), organism = 'mmusculus', 
                             ordered_query = T, evcodes = TRUE)
go_cm_GR <- gprofiler2::gost(query = list(Upregulated = cm_GR_up$Gene_names,
                                          Downregulated = cm_GR_down$Gene_names), 
                             sources = c("GO:BP", "GO:CC", "GO:MF"), organism = 'mmusculus', 
                             ordered_query = T, evcodes = TRUE)
#View(go_cm_all$result)
#View(go_cm_PR_up$result)
#View(go_cm_PR_down$result)

go_cm_all$result %>%
  #dplyr::filter(-log10(p_value) > 15) %>%
  arrange(source, desc(query), desc(p_value)) %>%
  mutate(term_code = paste0(term_name, " (", term_id, ")")) %>%
  mutate(term_code=factor(term_code, levels=unique(term_code))) %>%
  ggplot(aes(y = term_code, group = query, x = -log10(p_value), fill = query)) +
  geom_col(position = "dodge", width = 0.8) +
  #geom_point(size = 10, show.legend = F) +
  geom_text(aes(label = intersection_size), 
            nudge_x = -1, size = 0.8) +
 # facet_grid(facets = vars(source), 
#             scales = "free_y", space = "free_y") +
  theme_classic() +
  scale_fill_manual(values = c("#FAD510", "#3B9AB2")) +
  labs(title = "Enriched GO terms for significant genes in both mouse models", 
       y = "Term name and code", 
       x = expression(paste("-lo", g[10], " (", italic("P"), " value)")), fill = "") + 
  theme(legend.position = "top",
        legend.justification='left',
        legend.key.size = unit(0.5,"line"),
        strip.text = element_text(
          size = 4)) +
  theme(text = element_text(size = 4)) +
  theme(plot.title = element_text(hjust = -0.5, size = 8)) +
  theme(legend.text=element_text(size=4))
#ggsave("~/Desktop/GO_both.pdf", width = 5, height = 4.25)

go_cm_PR$result %>%
  #dplyr::filter(-log10(p_value) > 20) %>%
  arrange(source, desc(query), desc(p_value)) %>%
  mutate(term_code = paste0(term_name, " (", term_id, ")")) %>%
  mutate(term_code=factor(term_code, levels=unique(term_code))) %>%
  ggplot(aes(y = term_code, group = query, x = -log10(p_value), fill = query)) +
  geom_col(position = "dodge", width = 0.8) +
  #geom_point(size = 10, show.legend = F) +
  geom_text(aes(label = intersection_size), 
            nudge_x = -1, size = 0.8) +
 # facet_grid(facets = vars(source), 
#             scales = "free_y", space = "free_y") +
  theme_classic() +
  scale_fill_manual(values = c("#FAD510", "#3B9AB2")) +
  labs(title = "Enriched GO terms for significant genes in (PR)400 mouse model", 
       y = "Term name and code", 
       x = expression(paste("-lo", g[10], " (", italic("P"), " value)")), fill = "") + 
  theme(legend.position = "top",
        legend.justification='left',
        legend.key.size = unit(0.5,"line"),
        strip.text = element_text(
          size = 4)) +
  theme(text = element_text(size = 4)) +
  theme(plot.title = element_text(hjust = -0.5, size = 8)) +
  theme(legend.text=element_text(size=4))
#ggsave("~/Desktop/GO_PR.pdf", width = 5, height = 4.25)

go_cm_GR$result %>%
  #dplyr::filter(-log10(p_value) > 20) %>%
  arrange(source, desc(query), desc(p_value)) %>%
  mutate(term_code = paste0(term_name, " (", term_id, ")")) %>%
  mutate(term_code=factor(term_code, levels=unique(term_code))) %>%
  ggplot(aes(y = term_code, group = query, x = -log10(p_value), fill = query)) +
  geom_col(position = "dodge", width = 0.8) +
  #geom_point(size = 10, show.legend = F) +
  geom_text(aes(label = intersection_size), 
            nudge_x = -1, size = 0.8) +
#  facet_grid(facets = vars(source), 
#             scales = "free_y", space = "free_y") +
  theme_classic() +
  scale_fill_manual(values = c("#FAD510", "#3B9AB2")) +
  labs(title = "Enriched GO terms for significant genes in (GR)400 mouse model", 
       y = "Term name and code", 
       x = expression(paste("-lo", g[10], " (", italic("P"), " value)")), fill = "") + 
  theme(legend.position = "top",
        legend.justification='left',
        legend.key.size = unit(0.5,"line"),
        strip.text = element_text(
          size = 4)) +
  theme(text = element_text(size = 4)) +
  theme(plot.title = element_text(hjust = -0.5, size = 8)) +
  theme(legend.text=element_text(size=4))
#ggsave("~/Desktop/GO_GR.pdf", width = 5, height = 4.25)


###venn diagram - check numbers first
up_PR_only <- carmm %>% dplyr::filter(FC_PR_WT > 1 & p_value_PR_WT < 0.05 &
                                        (FC_GR_WT < 1 | p_value_GR_WT > 0.05)) 
up_GR_only <- carmm %>% dplyr::filter(FC_GR_WT > 1 & p_value_GR_WT < 0.05 &
                                        (FC_PR_WT < 1 | p_value_PR_WT > 0.05)) 
down_PR_only <- carmm %>% dplyr::filter(FC_PR_WT < 1 & p_value_PR_WT < 0.05 &
                                          (FC_GR_WT > 1 | p_value_GR_WT > 0.05))
down_GR_only <- carmm %>% dplyr::filter(FC_GR_WT < 1 & p_value_GR_WT < 0.05 &
                                          (FC_PR_WT > 1 | p_value_PR_WT > 0.05))
#PR_down 509 - all_down 250 - PR_only 258
#GR_down 474 - all_down 250 - GR_only 223
#PR_up 329 - all_up 214 - PR_only 115
#GR_up 309 - all_up 214 - GR_only 95

mats_down <- c("(GR)400" = 223, "(PR)400" = 258, "(PR)400&(GR)400" = 250)

fit_down <- euler(mats_down)
#pdf("Venn_down.pdf", width = 3.5, height = 2, pointsize = 12)
plot(fit_down, main = "Venn diagram showing downregulated genes in the two models",
     quantities = list(type = c("counts")), labels = FALSE,
     fills = c("#ffa040", "#55a0fb"))
#dev.off()

mats_up <- c("(GR)400" = 95, "(PR)400" = 115, "(PR)400&(GR)400" = 214)
fit_up <- euler(mats_up)
#pdf("Venn_up.pdf", width = 3.5, height = 2, pointsize = 12)
plot(fit_up, main = "Venn diagram showing upregulated genes in the two models",
     quantities = list(type = c("counts")), labels = FALSE,
     fills = c("#ffa040", "#55a0fb"))
#dev.off()



#####create list of genes according to GO category
View(go_cm_all$result)
go_cm_all$result <- go_cm_all$result[order(go_cm_all$result$p_value),]

###synapse GO:0045202
list_synaptic1 <- go_cm_all$result %>%
  dplyr::filter(term_id == "GO:0045202") %>% ####arbitrary
  dplyr::select(intersection) %>%
  pull()
list_synaptic <- unique(unlist(strsplit(paste(list_synaptic1, collapse = ","), ",")))
df_synaptic <- carmm %>% dplyr::filter(Gene_names %in% list_synaptic)
  #left_join(go_cm_PR_up$result[10,], )
write.table(df_synaptic, "synapse.tsv", sep = "\t", quote = F, row.names = F)

###mitochondrion GO:0005739
list_mito1 <- go_cm_all$result %>%
  dplyr::filter(term_id == "GO:0005739") %>% ####arbitrary
  dplyr::select(intersection) %>%
  pull()
list_mito <- unique(unlist(strsplit(paste(list_mito1, collapse = ","), ",")))
df_mito <- carmm %>% dplyr::filter(Gene_names %in% list_mito)
write.table(df_mito, "mito.tsv", sep = "\t", quote = F, row.names = F)

###ECM GO:0031012
list_ecm1 <- go_cm_all$result %>%
  dplyr::filter(term_id == "GO:0031012") %>% ####arbitrary
  dplyr::select(intersection) %>%
  pull()
list_ecm <- unique(unlist(strsplit(paste(list_ecm1, collapse = ","), ",")))
df_ecm <- carmm %>% dplyr::filter(Gene_names %in% list_ecm)
write.table(df_ecm, "ecm.tsv", sep = "\t", quote = F, row.names = F)



#carmm %>% dplyr::filter(Gene_names %in% c("Lama2","Lama4","Lamb2","Utrn","P2rx7","Myh9")) -> NMJ
