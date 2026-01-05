#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(optparse)
    library(jsonlite)
})

option_list <- list(
    make_option(c("--matrix"), type = "character"),
    make_option(c("--clusters"), type = "character"),
    make_option(c("--samplesheet"), type = "character"),
    make_option(c("--treated"), type = "character"),
    make_option(c("--control"), type = "character"),
    make_option(c("--group"), type = "character", default = "NA"),
    make_option(c("--fdr"), type = "double", default = 0.05),
    make_option(c("--lfc"), type = "double", default = 1.0),
    make_option(c("--bootstrap"), type = "integer", default = 1000),
    make_option(c("--k_value"), type = "integer", default = 100000)
)

opt <- parse_args(OptionParser(option_list = option_list))

matrix_df <- read.table(opt$matrix, header = TRUE, sep = "\t", check.names = FALSE)
clusters_df <- read.table(opt$clusters, header = TRUE, sep = "\t", check.names = FALSE)
samples_df <- read.csv(opt$samplesheet, stringsAsFactors = FALSE)

sample_ids <- samples_df$sample_id
conditions <- samples_df$condition

if (!all(sample_ids %in% colnames(matrix_df))) {
    stop("Matrix columns do not include all sample IDs")
}

counts <- as.matrix(matrix_df[, sample_ids, drop = FALSE])
rownames(counts) <- paste(matrix_df$chr, matrix_df$start, matrix_df$end, sep = ":")

# Differential testing
log2fc <- log2(rowMeans(counts[, conditions == opt$treated, drop = FALSE]) /
               rowMeans(counts[, conditions == opt$control, drop = FALSE]))

use_rots <- TRUE
tryCatch({
    suppressPackageStartupMessages(library(ROTS))
}, error = function(e) {
    use_rots <<- FALSE
})

if (use_rots) {
    groups <- ifelse(conditions == opt$treated, 1, 0)
    rots_res <- ROTS::ROTS(counts, groups = groups, B = opt$bootstrap, K = opt$k_value)
    pvals <- rots_res$pvalue
} else {
    pvals <- apply(counts, 1, function(x) {
        t.test(x[conditions == opt$treated], x[conditions == opt$control])$p.value
    })
}

fdr <- p.adjust(pvals, method = "BH")

results <- data.frame(
    chr = matrix_df$chr,
    start = matrix_df$start,
    end = matrix_df$end,
    log2FC = log2fc,
    pvalue = pvals,
    FDR = fdr
)

if (all(c("chr", "start", "end", "cluster") %in% colnames(clusters_df))) {
    clusters_key <- paste(clusters_df$chr, clusters_df$start, clusters_df$end, sep = ":")
    results$cluster <- clusters_df$cluster[match(paste(results$chr, results$start, results$end, sep = ":"), clusters_key)]
}

write.table(results, "chipbinner.differential.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

sig <- results[results$FDR <= opt$fdr & abs(results$log2FC) >= opt$lfc, , drop = FALSE]
up <- sig[sig$log2FC > 0, , drop = FALSE]
down <- sig[sig$log2FC < 0, , drop = FALSE]

write.table(up[, c("chr", "start", "end")], "chipbinner.significant_up.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(down[, c("chr", "start", "end")], "chipbinner.significant_down.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

summary <- data.frame(
    method = "chipbinner",
    group = opt$group,
    caller = "NA",
    treated_condition = opt$treated,
    control_condition = opt$control,
    n_tested = nrow(results),
    n_fdr_pass = nrow(sig),
    n_up = nrow(up),
    n_down = nrow(down)
)
write.table(summary, "chipbinner.summary.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# Plots
counts_log <- log2(counts + 1)
try({
    pca <- prcomp(t(counts_log), scale. = TRUE)
    pdf("plots/PCA.pdf")
    plot(pca$x[,1], pca$x[,2], col = as.factor(conditions), pch = 19,
         xlab = "PC1", ylab = "PC2", main = "PCA")
    legend("topright", legend = levels(as.factor(conditions)), col = 1:length(unique(conditions)), pch = 19)
    dev.off()
}, silent = TRUE)

try({
    pdf("plots/correlation_heatmap.pdf")
    heatmap(cor(counts_log), main = "Correlation")
    dev.off()
}, silent = TRUE)

try({
    treated_means <- rowMeans(counts[, conditions == opt$treated, drop = FALSE])
    control_means <- rowMeans(counts[, conditions == opt$control, drop = FALSE])
    pdf("plots/density_scatter.pdf")
    plot(control_means, treated_means, pch = 16, cex = 0.5,
         xlab = "Control mean", ylab = "Treated mean", main = "Density scatter")
    dev.off()
}, silent = TRUE)
