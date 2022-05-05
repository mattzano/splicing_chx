#1. First perform low count filter to remove isoforms with low expressiom
# Just pick isoforms that have an library size adjusted mean count of at least n reads in either of the conditions
# Using DESeq2's normalisation method (so count = count / size factor for sample)
counts <- featureCounts[,c(1,6:9,14:17)] %>%
  tibble::column_to_rownames('ensgene')
min_mean_reads <- 5
sample_tbl <- data.frame(sample_name = colnames(counts),
                   cond_col = rep(c("base_cond", "treat_cond"), each = 4))

message(glue::glue("Filtering isoforms for minimum mean count in any condition of at least - {min_mean_reads}"))

# Named vector of size factors for each sample in matrix
sfs <- DESeq2::estimateSizeFactorsForMatrix(counts)


# Divide each count column (sample) by its corresponding size factor
norm_counts <- sweep(x = counts,
                     MARGIN = 2, # operate on columns
                     STATS = sfs,
                     FUN = '/')

# Get list of samples (column names) for each condition
conds <- list(sample_tbl[sample_tbl[, 2] == "base_cond", "sample_name"],
              sample_tbl[sample_tbl[, 2] == "treat_cond", "sample_name"]
)


names(conds) <- c(base_cond, treat_cond)


# Calculate means for each condition      
# Matrix with cols base_cond | treat_cond (mean count)
means_norm_counts <- sapply(X = conds, FUN = function(x) {rowMeans(norm_counts[, x])})

# Get the max condition mean count for each sample
# Just a vector of max counts
maxs <- rowMaxs(means_norm_counts)

# Filter original count matrices for rows with a condition mean > 10
counts <- counts[maxs > min_mean_reads, ]
