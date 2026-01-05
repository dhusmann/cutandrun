#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(optparse)
})

option_list <- list(
    make_option(c("--samplesheet"), type = "character"),
    make_option(c("--treated"), type = "character"),
    make_option(c("--control"), type = "character"),
    make_option(c("--group"), type = "character", default = "NA"),
    make_option(c("--caller"), type = "character", default = "NA"),
    make_option(c("--fdr"), type = "double", default = 0.05),
    make_option(c("--lfc"), type = "double", default = 1.0),
    make_option(c("--min_overlap"), type = "integer", default = 2),
    make_option(c("--summits"), type = "integer", default = 0),
    make_option(c("--backend"), type = "character", default = "DESeq2"),
    make_option(c("--use_spikein"), type = "character", default = "false"),
    make_option(c("--prefix"), type = "character", default = "diffbind"),
    make_option(c("--extra_params"), type = "character", default = NULL)
)

opt <- parse_args(OptionParser(option_list = option_list))

suppressPackageStartupMessages({
    library(DiffBind)
})

use_spikein <- tolower(opt$use_spikein) %in% c("true", "1", "yes")

samplesheet <- read.csv(opt$samplesheet, stringsAsFactors = FALSE)
required_cols <- c("SampleID", "Condition", "bamReads", "Peaks", "Replicate")
missing <- setdiff(required_cols, colnames(samplesheet))
if (length(missing) > 0) {
    stop("Missing required columns in samplesheet: ", paste(missing, collapse = ", "))
}

samplesheet <- samplesheet[samplesheet$Condition %in% c(opt$treated, opt$control), , drop = FALSE]
if (nrow(samplesheet) == 0) {
    stop("No samples remain after filtering to contrast conditions")
}

sample_ids <- samplesheet$SampleID
scale_factors <- NULL
if ("SpikeinScaleFactor" %in% colnames(samplesheet)) {
    scale_factors <- suppressWarnings(as.numeric(samplesheet$SpikeinScaleFactor))
}

if (use_spikein && (is.null(scale_factors) || any(is.na(scale_factors)))) {
    message("Spike-in scaling requested but scale factors missing; falling back to default normalization")
    use_spikein <- FALSE
}

extra <- list()
if (!is.null(opt$extra_params) && file.exists(opt$extra_params)) {
    suppressPackageStartupMessages({
        library(jsonlite)
        library(yaml)
    })
    if (grepl("\\.ya?ml$", opt$extra_params, ignore.case = TRUE)) {
        extra <- yaml::read_yaml(opt$extra_params)
    } else {
        extra <- jsonlite::fromJSON(opt$extra_params)
    }
}

backend <- toupper(opt$backend)
backend_method <- if (backend == "EDGER") DBA_EDGER else DBA_DESEQ2

summits_value <- if (opt$summits > 0) opt$summits else FALSE

count_args <- list(summits = summits_value, minOverlap = opt$min_overlap)
if (!is.null(extra$dba_count) && is.list(extra$dba_count)) {
    count_args <- modifyList(count_args, extra$dba_count)
}

analyze_args <- list(method = backend_method)
if (!is.null(extra$dba_analyze) && is.list(extra$dba_analyze)) {
    analyze_args <- modifyList(analyze_args, extra$dba_analyze)
}

contrast_args <- list(categories = DBA_CONDITION, group1 = opt$treated, group2 = opt$control)
if (!is.null(extra$dba_contrast) && is.list(extra$dba_contrast)) {
    contrast_args <- modifyList(contrast_args, extra$dba_contrast)
}

extract_counts <- function(dba_obj, sample_ids) {
    df <- NULL
    counts <- NULL
    peaks_df <- NULL
    try({
        df <- dba.peakset(dba_obj, bRetrieve = TRUE, DataType = DBA_DATA_FRAME)
    }, silent = TRUE)
    if (!is.null(df)) {
        count_cols <- match(sample_ids, colnames(df))
        if (all(!is.na(count_cols))) {
            counts <- as.matrix(df[, count_cols, drop = FALSE])
        }
        coord_cols <- list(chr = NULL, start = NULL, end = NULL)
        for (col in colnames(df)) {
            if (tolower(col) %in% c("chr", "chrom", "seqnames", "seqname")) coord_cols$chr <- col
            if (tolower(col) %in% c("start")) coord_cols$start <- col
            if (tolower(col) %in% c("end")) coord_cols$end <- col
        }
        if (!is.null(coord_cols$chr) && !is.null(coord_cols$start) && !is.null(coord_cols$end)) {
            peaks_df <- df[, c(coord_cols$chr, coord_cols$start, coord_cols$end), drop = FALSE]
            colnames(peaks_df) <- c("chr", "start", "end")
        }
    }
    if (is.null(counts)) {
        try({
            gr <- dba.peakset(dba_obj, bRetrieve = TRUE)
            df_gr <- as.data.frame(gr)
            count_cols <- match(sample_ids, colnames(df_gr))
            if (all(!is.na(count_cols))) {
                counts <- as.matrix(df_gr[, count_cols, drop = FALSE])
            }
            if (all(c("seqnames", "start", "end") %in% colnames(df_gr))) {
                peaks_df <- df_gr[, c("seqnames", "start", "end")]
                colnames(peaks_df) <- c("chr", "start", "end")
            }
        }, silent = TRUE)
    }
    list(counts = counts, peaks = peaks_df)
}

run_deseq2 <- function(counts, conditions, scale_factors, treated, control) {
    suppressPackageStartupMessages(library(DESeq2))
    cond <- factor(conditions, levels = c(control, treated))
    dds <- DESeqDataSetFromMatrix(countData = counts, colData = data.frame(condition = cond), design = ~condition)
    if (!is.null(scale_factors)) {
        sizeFactors(dds) <- 1 / scale_factors
    }
    dds <- DESeq(dds)
    res <- results(dds, contrast = c("condition", treated, control))
    data.frame(log2FC = res$log2FoldChange, pvalue = res$pvalue, FDR = res$padj)
}

run_edger <- function(counts, conditions, scale_factors, treated, control) {
    suppressPackageStartupMessages(library(edgeR))
    cond <- factor(conditions, levels = c(control, treated))
    y <- DGEList(counts = counts, group = cond)
    if (!is.null(scale_factors)) {
        y$samples$norm.factors <- 1 / scale_factors
    } else {
        y <- calcNormFactors(y)
    }
    design <- model.matrix(~cond)
    y <- estimateDisp(y, design)
    fit <- glmQLFit(y, design)
    qlf <- glmQLFTest(fit, coef = 2)
    data.frame(log2FC = qlf$table$logFC, pvalue = qlf$table$PValue, FDR = p.adjust(qlf$table$PValue, method = "BH"))
}

write_plot_placeholder <- function(path, title) {
    pdf(path)
    plot.new()
    title(main = title)
    dev.off()
}

samplesheet$Condition <- as.character(samplesheet$Condition)

if (use_spikein) {
    dba_obj <- dba(sampleSheet = samplesheet)
    dba_obj <- do.call(dba.count, c(list(dba_obj), count_args))
    dba_obj <- do.call(dba.contrast, c(list(dba_obj), contrast_args))

    extracted <- extract_counts(dba_obj, sample_ids)
    if (is.null(extracted$counts) || is.null(extracted$peaks)) {
        message("Unable to extract counts/peaks for spike-in normalization; falling back to DiffBind default")
        use_spikein <- FALSE
    } else {
        counts <- extracted$counts
        if (backend == "EDGER") {
            res_stats <- run_edger(counts, samplesheet$Condition, scale_factors, opt$treated, opt$control)
        } else {
            res_stats <- run_deseq2(counts, samplesheet$Condition, scale_factors, opt$treated, opt$control)
        }
        res_df <- cbind(extracted$peaks, res_stats)
        res_df$log2FC[is.na(res_df$log2FC)] <- 0
        res_df$FDR[is.na(res_df$FDR)] <- 1
        res_df$pvalue[is.na(res_df$pvalue)] <- 1

        write.table(res_df, paste0(opt$prefix, ".results.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

        sig <- res_df[res_df$FDR <= opt$fdr & abs(res_df$log2FC) >= opt$lfc, , drop = FALSE]
        write.table(sig[, c("chr", "start", "end")], paste0(opt$prefix, ".significant.bed"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
        up <- sig[sig$log2FC > 0, , drop = FALSE]
        down <- sig[sig$log2FC < 0, , drop = FALSE]
        write.table(up[, c("chr", "start", "end")], paste0(opt$prefix, ".significant_up.bed"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
        write.table(down[, c("chr", "start", "end")], paste0(opt$prefix, ".significant_down.bed"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

        n_tested <- nrow(res_df)
        n_fdr <- nrow(sig)
        n_up <- nrow(up)
        n_down <- nrow(down)
        summary <- data.frame(
            method = "diffbind",
            group = opt$group,
            caller = opt$caller,
            treated_condition = opt$treated,
            control_condition = opt$control,
            n_tested = n_tested,
            n_fdr_pass = n_fdr,
            n_up = n_up,
            n_down = n_down
        )
        write.table(summary, paste0(opt$prefix, ".summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

        # plots
        counts_log <- log2(counts + 1)
        try({
            pca <- prcomp(t(counts_log), scale. = TRUE)
            pdf("plots/PCA.pdf")
            plot(pca$x[,1], pca$x[,2], col = as.factor(samplesheet$Condition), pch = 19,
                 xlab = "PC1", ylab = "PC2", main = "PCA")
            legend("topright", legend = levels(as.factor(samplesheet$Condition)), col = 1:length(unique(samplesheet$Condition)), pch = 19)
            dev.off()
        }, silent = TRUE)
        try({
            pdf("plots/correlation_heatmap.pdf")
            heatmap(cor(counts_log), main = "Correlation")
            dev.off()
        }, silent = TRUE)
        try({
            pdf("plots/MA.pdf")
            plot(rowMeans(counts_log), res_df$log2FC, pch = 16, cex = 0.5, main = "MA", xlab = "Mean", ylab = "log2FC")
            dev.off()
        }, silent = TRUE)
        try({
            pdf("plots/volcano.pdf")
            plot(res_df$log2FC, -log10(res_df$FDR), pch = 16, cex = 0.5, main = "Volcano", xlab = "log2FC", ylab = "-log10(FDR)")
            dev.off()
        }, silent = TRUE)
        try({
            pdf("plots/heatmap.pdf")
            heatmap(counts_log, Rowv = NA, Colv = NA, scale = "row", main = "Binding heatmap")
            dev.off()
        }, silent = TRUE)
    }
}

if (!use_spikein) {
    dba_obj <- dba(sampleSheet = samplesheet)
    dba_obj <- do.call(dba.count, c(list(dba_obj), count_args))
    dba_obj <- do.call(dba.contrast, c(list(dba_obj), contrast_args))
    dba_obj <- do.call(dba.analyze, c(list(dba_obj), analyze_args))

    res <- dba.report(dba_obj, th = 1, fold = 0)
    res_df <- as.data.frame(res)

    coord_cols <- list(chr = NULL, start = NULL, end = NULL)
    for (col in colnames(res_df)) {
        if (tolower(col) %in% c("chr", "chrom", "seqnames", "seqname")) coord_cols$chr <- col
        if (tolower(col) %in% c("start")) coord_cols$start <- col
        if (tolower(col) %in% c("end")) coord_cols$end <- col
    }
    if (!is.null(coord_cols$chr)) {
        res_df$chr <- res_df[[coord_cols$chr]]
    }
    if (!is.null(coord_cols$start)) {
        res_df$start <- res_df[[coord_cols$start]]
    }
    if (!is.null(coord_cols$end)) {
        res_df$end <- res_df[[coord_cols$end]]
    }
    if (!"chr" %in% colnames(res_df)) {
        res_df$chr <- NA
    }
    if (!"start" %in% colnames(res_df)) {
        res_df$start <- NA
    }
    if (!"end" %in% colnames(res_df)) {
        res_df$end <- NA
    }

    if (!"log2FC" %in% colnames(res_df)) {
        if ("Fold" %in% colnames(res_df)) {
            res_df$log2FC <- res_df$Fold
        } else if ("fold" %in% colnames(res_df)) {
            res_df$log2FC <- res_df$fold
        } else {
            res_df$log2FC <- NA_real_
        }
    }
    if (!"FDR" %in% colnames(res_df)) {
        res_df$FDR <- NA_real_
    }

    write.table(res_df, paste0(opt$prefix, ".results.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

    sig <- res_df[!is.na(res_df$FDR) & res_df$FDR <= opt$fdr & abs(res_df$log2FC) >= opt$lfc, , drop = FALSE]
    write.table(sig[, c("chr", "start", "end")], paste0(opt$prefix, ".significant.bed"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
    up <- sig[sig$log2FC > 0, , drop = FALSE]
    down <- sig[sig$log2FC < 0, , drop = FALSE]
    write.table(up[, c("chr", "start", "end")], paste0(opt$prefix, ".significant_up.bed"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
    write.table(down[, c("chr", "start", "end")], paste0(opt$prefix, ".significant_down.bed"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

    n_tested <- nrow(res_df)
    n_fdr <- nrow(sig)
    n_up <- nrow(up)
    n_down <- nrow(down)
    summary <- data.frame(
        method = "diffbind",
        group = opt$group,
        caller = opt$caller,
        treated_condition = opt$treated,
        control_condition = opt$control,
        n_tested = n_tested,
        n_fdr_pass = n_fdr,
        n_up = n_up,
        n_down = n_down
    )
    write.table(summary, paste0(opt$prefix, ".summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

    # plots
    try({
        pdf("plots/PCA.pdf")
        dba.plotPCA(dba_obj, DBA_CONDITION, label = DBA_ID)
        dev.off()
    }, silent = TRUE)

    try({
        pdf("plots/correlation_heatmap.pdf")
        dba.plotHeatmap(dba_obj, correlations = TRUE)
        dev.off()
    }, silent = TRUE)

    try({
        pdf("plots/MA.pdf")
        dba.plotMA(dba_obj)
        dev.off()
    }, silent = TRUE)

    try({
        pdf("plots/volcano.pdf")
        dba.plotVolcano(dba_obj)
        dev.off()
    }, silent = TRUE)

    try({
        pdf("plots/heatmap.pdf")
        dba.plotHeatmap(dba_obj)
        dev.off()
    }, silent = TRUE)
}

if (!file.exists(paste0(opt$prefix, ".results.tsv"))) {
    stop("DiffBind failed to produce results")
}

plots <- c("plots/PCA.pdf", "plots/correlation_heatmap.pdf", "plots/MA.pdf", "plots/volcano.pdf", "plots/heatmap.pdf")
for (plot_file in plots) {
    if (!file.exists(plot_file)) {
        write_plot_placeholder(plot_file, basename(plot_file))
    }
}
