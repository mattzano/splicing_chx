full_counts <- featureCounts[,1:9]
full_counts <- column_to_rownames(full_counts, "Geneid")
pca <- prcomp(t(full_counts),scale=F)
plot(pca$x[,1], pca$x[,2])
pca.var <- pca$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100,1)
barplot(pca.var.per, main="Scree Plot", xlab= "Principal Component", ylab= "Percent Variation")

library(ggplot2)
pca.data <- data.frame(Sample=rownames(pca$x), X=pca$x[,1], Y=pca$x[,2])
ggplot(data = pca.data, aes(x=X, y=Y, label=Sample)) +
  geom_point(alpha = 0.5) + geom_text() + 
  xlab(paste("PC1 - ", pca.var.per[1], "%", sep="")) +
  ylab(paste("PC2 - ", pca.var.per[2], "%", sep= "")) +
  theme_bw() +
  ggtitle("PCA Graph")


loading_scores <- pca$rotation[,1]
gene_scores <- abs(loading_scores)
gene_score_ranked <- sort(gene_scores, decreasing = TRUE)
top_100_genes <- names(gene_score_ranked[1:100])
pca$rotation[top_100_genes,1]
