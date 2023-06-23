# Volcano Plot Function

make_volcano_plot <- function(res) {
  res <- res %>% 
    as.data.frame() %>% 
    mutate(graph_alpha = ifelse((padj < 0.05) | (is.na(padj)), 1,0.2)) %>%
    mutate(graph_color = case_when(graph_alpha == 1 & log2FoldChange < 0 ~ "Downregulated",
                                   graph_alpha == 1 & log2FoldChange > 0 ~ "Upregulated",
                                   TRUE ~ "Not significant")) 
   
  my_plot <- ggplot(res, aes(x = log2FoldChange,
                             y= -log10(padj),
                             fill = graph_color)) +
    geom_point(aes(alpha = graph_alpha), pch = 21) +
    #theme(legend.position = 'top') +
    ggpubr::theme_pubr() +
    scale_fill_brewer(palette="Set2") +
    ylab(expression(paste("-Lo", g[10], "P-value"))) +
    xlab(expression(paste("Lo", g[2], "Fold Change"))) +
    theme(text = element_text(size = 24)) +
    theme(legend.text = element_text(size = 22)) +
    geom_hline(yintercept = -log10(0.05)) +
    geom_vline(xintercept = c(0), linetype = "dotted") +
    guides(alpha = "none", size = "none", fill = "none") +
    scale_x_continuous(limits = c(-8,8))
  
    #scale_y_continuous(limits = c(-20, 350), breaks = c(0, 50, 100, 150, 200, 250, 300, 350))
  
  #plot(my_plot)

  return(my_plot)
}
  
# here are some helper functions to do extra stuff to the volcano plot
#this one will add the gene names of those of interest to the plot
add_name_to_plot <-  function(res, gene_name){
  
  first_plot <- make_volcano_plot(res)
  names(first_plot$data) <- c("Geneid", "baseMean", "log2FoldChange", "lfcSE", 
                        "stat", "pvalue", "padj", "ensgene","symbol", "graph_alpha", "graph_color")
  
  
  with_name <- first_plot$data %>% 
    as.data.frame() %>%    
    #rownames_to_column('ensgene') %>% 
    #left_join(annotables::grch38, by = 'symbol') %>%
    dplyr::filter(symbol %in% gene_name) #%>%
    #filter(-log10(padj) > log10padj) %>% 
    #filter(abs(log2FoldChange) > log2FoldCut)
  
  print(dim(with_name))
  named_plot <- first_plot +
    geom_text_repel(data = with_name, aes(x = log2FoldChange,
                                          y= -log10(padj), label = symbol), 
                    point.padding = 0.3,
                    min.segment.length = 0,
                    nudge_y = -0.2,
                    box.padding  = 3,
                    max.overlaps = Inf,
                    size=4, show.legend = F)
  
  #named_plot <- first_plot + 
  #  geom_text(data = with_name, aes(label = symbol))
  
  #plot(named_plot)
  return(named_plot)
}

#this one will return the volcano plot with added gene name labels based on significance

label_significant <- function(res, log2FoldCut = 1, log10padj = 1){
  first_plot <- make_volcano_plot(res)
  names(first_plot$data) <- c("Geneid", "baseMean", "log2FoldChange", "lfcSE", 
                              "stat", "pvalue", "padj", "ensgene","symbol", "graph_alpha", "graph_color")
  
  with_name <-first_plot$data %>% 
    as.data.frame() %>% 
    #rownames_to_column("ensgene") %>% 
    #left_join(annotables::grch38)  %>% 
    dplyr::filter(-log10(padj) > log10padj) %>% 
    dplyr::filter(abs(log2FoldChange) > log2FoldCut) #%>%
#    #filter(symbol %in% list_gene)
#  
  print(dim(with_name))
  named_plot <- first_plot +
    geom_text_repel(data = with_name, aes(x = log2FoldChange,
                                          y= -log10(padj), label = symbol), 
                    point.padding = 0.3,
                    nudge_y = 0.2,
                    min.segment.length = 0,
                    box.padding  = 2,
                    max.overlaps = Inf,
                    size=4, show.legend = F)
  
  return(named_plot)
}
