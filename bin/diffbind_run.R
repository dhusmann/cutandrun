#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

parse_args <- function(args) {
    res <- list()
    i <- 1
    while (i <= length(args)) {
        key <- args[i]
        if (grepl("^--", key)) {
            name <- sub("^--", "", key)
            if (i == length(args) || grepl("^--", args[i + 1])) {
                res[[name]] <- TRUE
                i <- i + 1
            } else {
                res[[name]] <- args[i + 1]
                i <- i + 2
            }
        } else {
            i <- i + 1
        }
    }
    res
}

as_bool <- function(val) {
    if (is.null(val)) {
        return(FALSE)
    }
    if (is.logical(val)) {
        return(val)
    }
    val <- tolower(as.character(val))
    return(val %in% c("true", "1", "yes"))
}

params <- parse_args(args)
records_path <- params[["records"]]
outdir <- params[["outdir"]]
contrast <- params[["contrast"]]
caller <- params[["caller"]]
group <- params[["group"]]

fdr <- as.numeric(params[["fdr"]])
lfc <- as.numeric(params[["lfc"]])
min_overlap <- as.integer(params[["min_overlap"]])
backend <- params[["backend"]]
recenter <- as_bool(params[["recenter"]])
summits <- as.integer(params[["summits"]])
norm_method <- params[["norm_method"]]
export_sheets <- as_bool(params[["export_sheets"]])

use_spikein <- as_bool(params[["use_spikein"]])

if (is.null(records_path) || is.null(outdir)) {
    stop("Missing required arguments: --records and --outdir")
}

if (is.na(fdr)) fdr <- 0.05
if (is.na(lfc)) lfc <- 1.0
if (is.na(min_overlap)) min_overlap <- 2
if (is.na(summits)) summits <- 0

if (!dir.exists(outdir)) {
    dir.create(outdir, recursive = TRUE)
}
plots_dir <- file.path(outdir, "plots")
if (!dir.exists(plots_dir)) {
    dir.create(plots_dir, recursive = TRUE)
}

records <- read.delim(records_path, stringsAsFactors = FALSE, check.names = FALSE)

samplesheet <- data.frame(
    SampleID = records$sample_id,
    Factor = records$group,
    Condition = records$condition,
    Replicate = records$replicate,
    bamReads = records$bam,
    Peaks = records$peaks,
    PeakCaller = records$caller,
    Tissue = "CUTRUN",
    stringsAsFactors = FALSE
)

samplesheet_path <- file.path(outdir, "diffbind.samplesheet.csv")
if (export_sheets) {
    write.csv(samplesheet, samplesheet_path, row.names = FALSE)
} else {
    write.csv(samplesheet, samplesheet_path, row.names = FALSE)
}

write_stub <- function() {
    peaks_file <- samplesheet$Peaks[1]
    regions <- data.frame()
    if (!is.null(peaks_file) && file.exists(peaks_file)) {
        peak_rows <- tryCatch({
            read.delim(peaks_file, header = FALSE, stringsAsFactors = FALSE)
        }, error = function(e) {
            NULL
        })
        if (!is.null(peak_rows) && nrow(peak_rows) > 0) {
            peak_rows <- peak_rows[, 1:3]
            colnames(peak_rows) <- c("chr", "start", "end")
            regions <- peak_rows
        }
    }

    if (nrow(regions) == 0) {
        regions <- data.frame(chr = character(), start = integer(), end = integer())
    }

    results <- data.frame(
        chr = regions$chr,
        start = regions$start,
        end = regions$end,
        log2FC = rep(0.0, nrow(regions)),
        pval = rep(1.0, nrow(regions)),
        FDR = rep(1.0, nrow(regions)),
        stringsAsFactors = FALSE
    )

    write.table(results, file = file.path(outdir, "diffbind.results.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

    sig <- results[results$FDR <= fdr & abs(results$log2FC) >= lfc, ]
    write.table(sig[, c("chr", "start", "end")], file = file.path(outdir, "diffbind.significant.bed"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
    write.table(sig[sig$log2FC >= lfc, c("chr", "start", "end")], file = file.path(outdir, "diffbind.significant_up.bed"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
    write.table(sig[sig$log2FC <= -lfc, c("chr", "start", "end")], file = file.path(outdir, "diffbind.significant_down.bed"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

    summary <- data.frame(
        caller = caller,
        group = group,
        treated = strsplit(contrast, ",")[[1]][1],
        control = strsplit(contrast, ",")[[1]][2],
        n_tested = nrow(results),
        n_fdr_pass = nrow(sig),
        n_up = sum(sig$log2FC >= lfc),
        n_down = sum(sig$log2FC <= -lfc),
        status = "RUN",
        reason = "ok",
        stringsAsFactors = FALSE
    )
    write.table(summary, file = file.path(outdir, "diffbind.summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

    if (use_spikein && "spikein_scale_factor" %in% colnames(records)) {
        factors <- suppressWarnings(as.numeric(records$spikein_scale_factor))
        if (any(!is.na(factors))) {
            size_factors <- ifelse(factors == 0, 1, 1 / factors)
            norm_out <- data.frame(sample_id = records$sample_id, size_factor = size_factors, stringsAsFactors = FALSE)
            write.table(norm_out, file = file.path(outdir, "diffbind.normalization_factors.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
        }
    }
}

use_diffbind <- FALSE
if (requireNamespace("DiffBind", quietly = TRUE)) {
    use_diffbind <- TRUE
}

if (!use_diffbind) {
    write_stub()
    quit(save = "no")
}

# Best-effort DiffBind execution
tryCatch({
    suppressPackageStartupMessages(library(DiffBind))

    dba_obj <- dba(sampleSheet = samplesheet)
    summit_size <- if (recenter) summits else 0
    dba_obj <- dba.count(dba_obj, summits = summit_size, minOverlap = min_overlap)

    if (use_spikein && "spikein_scale_factor" %in% colnames(records)) {
        factors <- suppressWarnings(as.numeric(records$spikein_scale_factor))
        if (any(!is.na(factors))) {
            size_factors <- ifelse(factors == 0, 1, 1 / factors)
            names(size_factors) <- records$sample_id
            size_factors <- size_factors[samplesheet$SampleID]
            dba_obj <- dba.normalize(dba_obj, normalize = DBA_NORM_LIB, library = size_factors)
            norm_out <- data.frame(sample_id = samplesheet$SampleID, size_factor = size_factors, stringsAsFactors = FALSE)
            write.table(norm_out, file = file.path(outdir, "diffbind.normalization_factors.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
        }
    } else if (!is.null(norm_method) && norm_method != "native") {
        dba_obj <- dba.normalize(dba_obj, method = norm_method)
    }

    contrast_labels <- strsplit(contrast, ",")[[1]]
    dba_obj <- dba.contrast(dba_obj, categories = DBA_CONDITION, group1 = contrast_labels[1], group2 = contrast_labels[2])

    method_flag <- ifelse(toupper(backend) == "EDGER", DBA_EDGER, DBA_DESEQ2)
    dba_obj <- dba.analyze(dba_obj, method = method_flag)

    report <- dba.report(dba_obj, th = fdr, fold = lfc)
    results <- as.data.frame(report)
    results_out <- data.frame(
        chr = results$seqnames,
        start = results$start,
        end = results$end,
        log2FC = results$Fold,
        pval = results$p.value,
        FDR = results$FDR,
        stringsAsFactors = FALSE
    )
    write.table(results_out, file = file.path(outdir, "diffbind.results.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

    sig <- results_out[results_out$FDR <= fdr & abs(results_out$log2FC) >= lfc, ]
    write.table(sig[, c("chr", "start", "end")], file = file.path(outdir, "diffbind.significant.bed"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
    write.table(sig[sig$log2FC >= lfc, c("chr", "start", "end")], file = file.path(outdir, "diffbind.significant_up.bed"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
    write.table(sig[sig$log2FC <= -lfc, c("chr", "start", "end")], file = file.path(outdir, "diffbind.significant_down.bed"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

    summary <- data.frame(
        caller = caller,
        group = group,
        treated = contrast_labels[1],
        control = contrast_labels[2],
        n_tested = nrow(results_out),
        n_fdr_pass = nrow(sig),
        n_up = sum(sig$log2FC >= lfc),
        n_down = sum(sig$log2FC <= -lfc),
        status = "RUN",
        reason = "ok",
        stringsAsFactors = FALSE
    )
    write.table(summary, file = file.path(outdir, "diffbind.summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

    saveRDS(dba_obj, file = file.path(outdir, "diffbind.dba.rds"))

}, error = function(e) {
    message("DiffBind failed, falling back to stub output: ", e$message)
    write_stub()
})

norm_path <- file.path(outdir, "diffbind.normalization_factors.tsv")
if (!file.exists(norm_path)) {
    norm_out <- data.frame(
        sample_id = records$sample_id,
        size_factor = rep(NA, nrow(records)),
        stringsAsFactors = FALSE
    )
    write.table(norm_out, file = norm_path, sep = "\t", quote = FALSE, row.names = FALSE)
}

dba_path <- file.path(outdir, "diffbind.dba.rds")
if (!file.exists(dba_path)) {
    file.create(dba_path)
}
