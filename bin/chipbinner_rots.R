#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(optparse)
    library(jsonlite)
})

option_list <- list(
    make_option(c("--matrix"), type = "character"),
    make_option(c("--clusters"), type = "character"),
    make_option(c("--grid_summary"), type = "character"),
    make_option(c("--norm_info"), type = "character"),
    make_option(c("--samplesheet"), type = "character"),
    make_option(c("--treated"), type = "character"),
    make_option(c("--control"), type = "character"),
    make_option(c("--group"), type = "character", default = "NA"),
    make_option(c("--fdr"), type = "double", default = 0.05),
    make_option(c("--lfc"), type = "double", default = 1.0),
    make_option(c("--bootstrap"), type = "integer", default = 1000),
    make_option(c("--k_value"), type = "integer", default = 100000),
    make_option(c("--lola_run"), type = "character", default = "false")
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

cluster_ids <- rep(NA, nrow(results))
if (all(c("chr", "start", "end", "cluster") %in% colnames(clusters_df))) {
    clusters_key <- paste(clusters_df$chr, clusters_df$start, clusters_df$end, sep = ":")
    cluster_ids <- clusters_df$cluster[match(paste(results$chr, results$start, results$end, sep = ":"), clusters_key)]
}
results$cluster_id <- cluster_ids

cluster_labels <- rep("stable", length(cluster_ids))
if (!all(is.na(cluster_ids))) {
    treated_idx <- conditions == opt$treated
    control_idx <- conditions == opt$control
    for (cid in unique(cluster_ids)) {
        if (is.na(cid)) {
            next
        }
        if (cid == -1) {
            cluster_labels[cluster_ids == cid] <- "noise"
            next
        }
        cluster_rows <- cluster_ids == cid
        treated_mean <- mean(rowMeans(counts[cluster_rows, treated_idx, drop = FALSE]))
        control_mean <- mean(rowMeans(counts[cluster_rows, control_idx, drop = FALSE]))
        if (is.nan(treated_mean)) treated_mean <- 0
        if (is.nan(control_mean)) control_mean <- 0
        log2fc_cluster <- log2((treated_mean + 1e-6) / (control_mean + 1e-6))
        if (abs(log2fc_cluster) < opt$lfc) {
            label <- "stable"
        } else if (log2fc_cluster > 0) {
            label <- "treated_enriched"
        } else {
            label <- "control_enriched"
        }
        cluster_labels[cluster_ids == cid] <- label
    }
}
results$cluster_label <- cluster_labels

write.table(results, "chipbinner.differential.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

sig <- results[results$FDR <= opt$fdr & abs(results$log2FC) >= opt$lfc, , drop = FALSE]
up <- sig[sig$log2FC > 0, , drop = FALSE]
down <- sig[sig$log2FC < 0, , drop = FALSE]

write.table(up[, c("chr", "start", "end")], "chipbinner.significant_up.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(down[, c("chr", "start", "end")], "chipbinner.significant_down.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

write.table(results[results$cluster_label == "control_enriched", c("chr", "start", "end")],
            "chipbinner.control_enriched.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(results[results$cluster_label == "treated_enriched", c("chr", "start", "end")],
            "chipbinner.treated_enriched.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(results[results$cluster_label == "stable", c("chr", "start", "end")],
            "chipbinner.stable.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(results[results$cluster_label == "noise", c("chr", "start", "end")],
            "chipbinner.noise.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

selected_row <- NULL
if (!is.null(opt$grid_summary) && file.exists(opt$grid_summary)) {
    grid <- read.table(opt$grid_summary, header = TRUE, sep = "\t", check.names = FALSE)
    if ("selected" %in% colnames(grid)) {
        selected_row <- grid[grid$selected == TRUE | grid$selected == "True" | grid$selected == "true", , drop = FALSE]
        if (nrow(selected_row) > 0) {
            selected_row <- selected_row[1, , drop = FALSE]
        }
    }
}

norm_info <- list(use_input = FALSE, use_ms_scaling = FALSE)
if (!is.null(opt$norm_info) && file.exists(opt$norm_info)) {
    try({
        norm_info <- jsonlite::fromJSON(opt$norm_info)
    }, silent = TRUE)
}

use_spikein_flag <- ifelse(is.null(norm_info$use_spikein_scaling), "false",
                           ifelse(norm_info$use_spikein_scaling, "true", "false"))

n_clusters <- length(unique(cluster_ids[!is.na(cluster_ids) & cluster_ids != -1]))
min_cluster_size <- if (!is.null(selected_row) && "min_cluster_size" %in% colnames(selected_row)) selected_row$min_cluster_size else NA
min_samples <- if (!is.null(selected_row) && "min_samples" %in% colnames(selected_row)) selected_row$min_samples else NA

summary <- data.frame(
    method = "chipbinner",
    group = opt$group,
    caller = "NA",
    treated_condition = opt$treated,
    control_condition = opt$control,
    n_bins_tested = nrow(results),
    n_tested = nrow(results),
    n_fdr_pass = nrow(sig),
    n_up = nrow(up),
    n_down = nrow(down),
    n_clusters = n_clusters,
    min_cluster_size = min_cluster_size,
    min_samples = min_samples,
    input_normalization_used = ifelse(is.null(norm_info$use_input), FALSE, norm_info$use_input),
    ms_scaling_used = ifelse(is.null(norm_info$use_ms_scaling), FALSE, norm_info$use_ms_scaling),
    use_spikein = use_spikein_flag,
    span_mode_used = "NA",
    lola_run = tolower(opt$lola_run) %in% c("true", "1", "yes")
)
write.table(summary, "chipbinner.summary.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# Plots
# counts already include pseudocount from normalization
counts_log <- log2(counts)
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
    cluster_factor <- as.factor(cluster_labels)
    cols <- as.numeric(cluster_factor)
    plot(control_means, treated_means, pch = 16, cex = 0.5, col = cols,
         xlab = "Control mean", ylab = "Treated mean", main = "Density scatter")
    legend("topleft", legend = levels(cluster_factor), col = seq_along(levels(cluster_factor)), pch = 16, cex = 0.6)
    dev.off()
}, silent = TRUE)
