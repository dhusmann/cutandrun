#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(optparse)
})

option_list <- list(
    make_option(c("--samplesheet"), type = "character"),
    make_option(c("--treated"), type = "character"),
    make_option(c("--control"), type = "character"),
    make_option(c("--group"), type = "character", default = "NA"),
    make_option(c("--fdr"), type = "double", default = 0.05),
    make_option(c("--backend"), type = "character", default = "DESeq2"),
    make_option(c("--use_spikein"), type = "character", default = "false"),
    make_option(c("--prefix"), type = "character", default = "span_fallback")
)

opt <- parse_args(OptionParser(option_list = option_list))

suppressPackageStartupMessages({
    library(GenomicRanges)
    library(GenomicAlignments)
    library(Rsamtools)
})

use_spikein <- tolower(opt$use_spikein) %in% c("true", "1", "yes")

samplesheet <- read.csv(opt$samplesheet, stringsAsFactors = FALSE)
required_cols <- c("SampleID", "Condition", "bamReads", "Peaks")
missing <- setdiff(required_cols, colnames(samplesheet))
if (length(missing) > 0) {
    stop("Missing required columns in samplesheet: ", paste(missing, collapse = ", "))
}

samplesheet <- samplesheet[samplesheet$Condition %in% c(opt$treated, opt$control), , drop = FALSE]
if (nrow(samplesheet) == 0) {
    stop("No samples remain after filtering to contrast conditions")
}

read_peaks <- function(path) {
    df <- tryCatch({
        read.table(path, sep = "\t", header = FALSE, comment.char = "", quote = "", stringsAsFactors = FALSE)
    }, error = function(e) {
        stop("Failed to read peaks file: ", path)
    })
    if (ncol(df) < 3) {
        stop("Peaks file has fewer than 3 columns: ", path)
    }
    df <- df[, 1:3]
    suppressWarnings({
        df <- df[!is.na(as.numeric(df[, 2])) & !is.na(as.numeric(df[, 3])), , drop = FALSE]
    })
    colnames(df) <- c("chr", "start", "end")
    df$start <- as.integer(df$start)
    df$end <- as.integer(df$end)
    df
}

if (any(is.na(samplesheet$Peaks)) || any(samplesheet$Peaks == "")) {
    stop("Missing SPAN peak files in samplesheet for fallback differential.")
}

peak_list <- lapply(samplesheet$Peaks, read_peaks)
peaks_list_path <- paste0(opt$prefix, ".peaks.list")
writeLines(paste(samplesheet$Peaks, collapse = ";"), peaks_list_path)
gr_list <- lapply(peak_list, function(df) {
    GRanges(seqnames = df$chr, ranges = IRanges(start = df$start + 1L, end = df$end))
})

union_gr <- reduce(do.call(c, gr_list))

union_df <- data.frame(
    chr = as.character(seqnames(union_gr)),
    start = start(union_gr) - 1L,
    end = end(union_gr),
    stringsAsFactors = FALSE
)

union_bed <- paste0(opt$prefix, ".union.bed")
write.table(union_df, union_bed, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

bam_files <- samplesheet$bamReads
sample_ids <- samplesheet$SampleID
names(bam_files) <- sample_ids

detect_paired <- function(bam) {
    bf <- BamFile(bam, yieldSize = 1000)
    open(bf)
    flags <- tryCatch({
        scanBam(bf, param = ScanBamParam(what = "flag"))[[1]]$flag
    }, error = function(e) {
        integer()
    })
    close(bf)
    if (length(flags) == 0) {
        return(FALSE)
    }
    any(bitwAnd(flags, 1L) != 0)
}

paired_flags <- vapply(bam_files, detect_paired, logical(1))
single_end <- !any(paired_flags)

bam_list <- BamFileList(bam_files, yieldSize = 1000000)
se <- summarizeOverlaps(
    features = union_gr,
    reads = bam_list,
    mode = "Union",
    singleEnd = single_end,
    ignore.strand = TRUE
)

counts <- assay(se)
colnames(counts) <- sample_ids

counts_df <- cbind(union_df, counts)
counts_path <- paste0(opt$prefix, ".counts.tsv")
write.table(counts_df, counts_path, sep = "\t", quote = FALSE, row.names = FALSE)

scale_factors <- NULL
if ("SpikeinScaleFactor" %in% colnames(samplesheet)) {
    scale_factors <- suppressWarnings(as.numeric(samplesheet$SpikeinScaleFactor))
}

if (use_spikein && (is.null(scale_factors) || any(is.na(scale_factors)))) {
    message("Spike-in scaling requested but scale factors missing; falling back to default normalization")
    use_spikein <- FALSE
}

backend <- toupper(opt$backend)

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

conditions <- samplesheet$Condition
if (backend == "EDGER") {
    res_stats <- run_edger(counts, conditions, if (use_spikein) scale_factors else NULL, opt$treated, opt$control)
} else {
    res_stats <- run_deseq2(counts, conditions, if (use_spikein) scale_factors else NULL, opt$treated, opt$control)
}

res_stats$log2FC[is.na(res_stats$log2FC)] <- 0
res_stats$pvalue[is.na(res_stats$pvalue)] <- 1
res_stats$FDR[is.na(res_stats$FDR)] <- 1

res_df <- cbind(union_df, res_stats)
write.table(res_df, paste0(opt$prefix, ".results.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

sig <- res_df[res_df$FDR <= opt$fdr, , drop = FALSE]
write.table(sig[, c("chr", "start", "end")], paste0(opt$prefix, ".significant.bed"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
up <- sig[sig$log2FC > 0, , drop = FALSE]
down <- sig[sig$log2FC < 0, , drop = FALSE]
write.table(up[, c("chr", "start", "end")], paste0(opt$prefix, ".significant_up.bed"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(down[, c("chr", "start", "end")], paste0(opt$prefix, ".significant_down.bed"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

summary <- data.frame(
    method = "span_fallback",
    group = opt$group,
    caller = "NA",
    treated_condition = opt$treated,
    control_condition = opt$control,
    n_tested = nrow(res_df),
    n_fdr_pass = nrow(sig),
    n_up = nrow(up),
    n_down = nrow(down),
    use_spikein = if (use_spikein) "true" else "false",
    span_mode_used = "fallback"
)
write.table(summary, paste0(opt$prefix, ".summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
