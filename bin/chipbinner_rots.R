#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(optparse)
    library(jsonlite)
})

option_list <- list(
    make_option(c("--matrix"), type = "character"),
    make_option(c("--clusters"), type = "character"),
    make_option(c("--clusters_2"), type = "character", default = ""),
    make_option(c("--clusters_3"), type = "character", default = ""),
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
samples_df <- read.csv(opt$samplesheet, stringsAsFactors = FALSE)

sample_ids <- samples_df$sample_id
conditions <- samples_df$condition
treated_idx <- conditions == opt$treated
control_idx <- conditions == opt$control

if (!all(sample_ids %in% colnames(matrix_df))) {
    stop("Matrix columns do not include all sample IDs")
}

counts <- as.matrix(matrix_df[, sample_ids, drop = FALSE])
rownames(counts) <- paste(matrix_df$chr, matrix_df$start, matrix_df$end, sep = ":")

# Differential testing
log2fc <- log2(rowMeans(counts[, treated_idx, drop = FALSE]) /
               rowMeans(counts[, control_idx, drop = FALSE]))

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
    if (sum(treated_idx) < 2 || sum(control_idx) < 2) {
        warning("ROTS unavailable and <2 samples per condition; assigning p-values of 1")
        pvals <- rep(1, nrow(counts))
    } else {
        pvals <- apply(counts, 1, function(x) {
            t.test(x[treated_idx], x[control_idx])$p.value
        })
    }
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

read_clusters <- function(path) {
    if (is.null(path) || path == "" || !file.exists(path)) {
        return(NULL)
    }
    df <- tryCatch(read.table(path, header = TRUE, sep = "\t", check.names = FALSE),
                   error = function(e) NULL)
    if (is.null(df) || nrow(df) == 0) {
        return(NULL)
    }
    if (!all(c("chr", "start", "end", "cluster") %in% colnames(df))) {
        return(NULL)
    }
    df
}

cluster_ids_from_df <- function(clusters_df, results_df) {
    if (is.null(clusters_df)) {
        return(NULL)
    }
    clusters_key <- paste(clusters_df$chr, clusters_df$start, clusters_df$end, sep = ":")
    cluster_ids <- clusters_df$cluster[match(paste(results_df$chr, results_df$start, results_df$end, sep = ":"), clusters_key)]
    cluster_ids
}

compute_cluster_labels <- function(cluster_ids, counts_mat, conds, treated, control, lfc) {
    if (is.null(cluster_ids) || length(cluster_ids) == 0) {
        return(NULL)
    }
    labels <- rep("stable", length(cluster_ids))
    if (all(is.na(cluster_ids))) {
        labels[] <- "NA"
        return(labels)
    }
    treated_idx <- conds == treated
    control_idx <- conds == control
    for (cid in unique(cluster_ids)) {
        if (is.na(cid)) {
            next
        }
        if (cid == -1) {
            labels[cluster_ids == cid] <- "noise"
            next
        }
        cluster_rows <- cluster_ids == cid
        treated_mean <- mean(rowMeans(counts_mat[cluster_rows, treated_idx, drop = FALSE]))
        control_mean <- mean(rowMeans(counts_mat[cluster_rows, control_idx, drop = FALSE]))
        if (is.nan(treated_mean)) treated_mean <- 0
        if (is.nan(control_mean)) control_mean <- 0
        log2fc_cluster <- log2((treated_mean + 1e-6) / (control_mean + 1e-6))
        if (abs(log2fc_cluster) < lfc) {
            label <- "stable"
        } else if (log2fc_cluster > 0) {
            label <- "treated_enriched"
        } else {
            label <- "control_enriched"
        }
        labels[cluster_ids == cid] <- label
    }
    labels
}

write_cluster_beds <- function(results_df, labels, suffix, model, manifest_df) {
    if (is.null(labels)) {
        return(manifest_df)
    }
    suffix_tag <- if (suffix == "") "" else paste0(".", suffix)
    out_files <- list(
        control_enriched = paste0("chipbinner.control_enriched", suffix_tag, ".bed"),
        treated_enriched = paste0("chipbinner.treated_enriched", suffix_tag, ".bed"),
        stable = paste0("chipbinner.stable", suffix_tag, ".bed"),
        noise = paste0("chipbinner.noise", suffix_tag, ".bed")
    )
    write.table(results_df[labels == "control_enriched", c("chr", "start", "end")],
                out_files$control_enriched, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
    write.table(results_df[labels == "treated_enriched", c("chr", "start", "end")],
                out_files$treated_enriched, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
    write.table(results_df[labels == "stable", c("chr", "start", "end")],
                out_files$stable, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
    write.table(results_df[labels == "noise", c("chr", "start", "end")],
                out_files$noise, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

    for (label in names(out_files)) {
        bed_file <- out_files[[label]]
        if (file.exists(bed_file)) {
            manifest_df <- rbind(
                manifest_df,
                data.frame(
                    model = model,
                    label = label,
                    bed = bed_file,
                    stringsAsFactors = FALSE
                )
            )
        }
    }
    manifest_df
}

clusters_best <- read_clusters(opt$clusters)
clusters_2 <- read_clusters(opt$clusters_2)
clusters_3 <- read_clusters(opt$clusters_3)

cluster_ids_best <- cluster_ids_from_df(clusters_best, results)
cluster_labels_best <- compute_cluster_labels(cluster_ids_best, counts, conditions, opt$treated, opt$control, opt$lfc)
if (is.null(cluster_ids_best)) {
    cluster_ids_best <- rep(NA, nrow(results))
}
if (is.null(cluster_labels_best)) {
    cluster_labels_best <- rep("NA", nrow(results))
}

cluster_ids_2_raw <- cluster_ids_from_df(clusters_2, results)
cluster_labels_2_raw <- compute_cluster_labels(cluster_ids_2_raw, counts, conditions, opt$treated, opt$control, opt$lfc)
cluster_ids_2 <- cluster_ids_2_raw
cluster_labels_2 <- cluster_labels_2_raw
if (is.null(cluster_ids_2)) {
    cluster_ids_2 <- rep(NA, nrow(results))
}
if (is.null(cluster_labels_2)) {
    cluster_labels_2 <- rep("NA", nrow(results))
}

cluster_ids_3_raw <- cluster_ids_from_df(clusters_3, results)
cluster_labels_3_raw <- compute_cluster_labels(cluster_ids_3_raw, counts, conditions, opt$treated, opt$control, opt$lfc)
cluster_ids_3 <- cluster_ids_3_raw
cluster_labels_3 <- cluster_labels_3_raw
if (is.null(cluster_ids_3)) {
    cluster_ids_3 <- rep(NA, nrow(results))
}
if (is.null(cluster_labels_3)) {
    cluster_labels_3 <- rep("NA", nrow(results))
}

results$cluster_id <- cluster_ids_best
results$cluster_label <- cluster_labels_best
results$cluster_id_2clusters <- cluster_ids_2
results$cluster_label_2clusters <- cluster_labels_2
results$cluster_id_3clusters <- cluster_ids_3
results$cluster_label_3clusters <- cluster_labels_3

write.table(results, "chipbinner.differential.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

sig <- results[results$FDR <= opt$fdr & abs(results$log2FC) >= opt$lfc, , drop = FALSE]
up <- sig[sig$log2FC > 0, , drop = FALSE]
down <- sig[sig$log2FC < 0, , drop = FALSE]

write.table(up[, c("chr", "start", "end")], "chipbinner.significant_up.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(down[, c("chr", "start", "end")], "chipbinner.significant_down.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

bed_manifest <- data.frame(
    model = character(),
    label = character(),
    bed = character(),
    stringsAsFactors = FALSE
)
bed_manifest <- write_cluster_beds(results, cluster_labels_best, "", "best", bed_manifest)

if (!is.null(cluster_labels_2_raw)) {
    bed_manifest <- write_cluster_beds(results, cluster_labels_2_raw, "2clusters", "2clusters", bed_manifest)
}

if (!is.null(cluster_labels_3_raw)) {
    bed_manifest <- write_cluster_beds(results, cluster_labels_3_raw, "3clusters", "3clusters", bed_manifest)
}

write.table(bed_manifest, "chipbinner.cluster_beds.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

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

n_clusters <- length(unique(cluster_ids_best[!is.na(cluster_ids_best) & cluster_ids_best != -1]))
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
    cluster_factor <- as.factor(cluster_labels_best)
    cols <- as.numeric(cluster_factor)
    plot(control_means, treated_means, pch = 16, cex = 0.5, col = cols,
         xlab = "Control mean", ylab = "Treated mean", main = "Density scatter")
    legend("topleft", legend = levels(cluster_factor), col = seq_along(levels(cluster_factor)), pch = 16, cex = 0.6)
    dev.off()
}, silent = TRUE)
