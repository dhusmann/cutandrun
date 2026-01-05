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

read_extra_params <- function(path) {
    if (is.null(path) || !nzchar(path)) {
        return(list())
    }
    if (!file.exists(path)) {
        stop(sprintf("extra_params file not found: %s", path))
    }
    ext <- tolower(tools::file_ext(path))
    if (ext %in% c("yml", "yaml")) {
        if (!requireNamespace("yaml", quietly = TRUE)) {
            stop("YAML extra_params requires the 'yaml' R package; install it or provide JSON instead.")
        }
        return(yaml::yaml.load_file(path))
    }
    if (requireNamespace("jsonlite", quietly = TRUE)) {
        return(jsonlite::fromJSON(path, simplifyVector = FALSE))
    }
    if (requireNamespace("rjson", quietly = TRUE)) {
        return(rjson::fromJSON(file = path))
    }
    stop("JSON extra_params requires the 'jsonlite' or 'rjson' R package.")
}

extra_params_for <- function(extra, key) {
    if (is.null(extra) || !is.list(extra)) {
        return(list())
    }
    if (!is.null(extra[[key]])) {
        return(extra[[key]])
    }
    key_dot <- gsub("_", ".", key)
    if (!is.null(extra[[key_dot]])) {
        return(extra[[key_dot]])
    }
    key_us <- gsub("\\.", "_", key)
    if (!is.null(extra[[key_us]])) {
        return(extra[[key_us]])
    }
    return(list())
}

resolve_constants <- function(value) {
    if (is.list(value)) {
        return(lapply(value, resolve_constants))
    }
    if (is.character(value) && length(value) == 1) {
        if (exists(value, envir = asNamespace("DiffBind"), inherits = FALSE)) {
            return(get(value, envir = asNamespace("DiffBind")))
        }
    }
    value
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
allow_partial <- as_bool(params[["allow_partial"]])
extra_params_path <- params[["extra_params"]]
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
write.csv(samplesheet, samplesheet_path, row.names = FALSE)

contrast_labels <- strsplit(contrast, ",")[[1]]
contrast_labels <- trimws(contrast_labels)
if (length(contrast_labels) != 2) {
    stop("contrast must contain exactly two comma-separated labels")
}

write_summary <- function(status, reason, n_tested = 0, n_fdr_pass = 0, n_up = 0, n_down = 0) {
    summary <- data.frame(
        caller = caller,
        group = group,
        treated = contrast_labels[1],
        control = contrast_labels[2],
        n_tested = n_tested,
        n_fdr_pass = n_fdr_pass,
        n_up = n_up,
        n_down = n_down,
        status = status,
        reason = reason,
        stringsAsFactors = FALSE
    )
    write.table(summary, file = file.path(outdir, "diffbind.summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
}

write_error <- function(message) {
    writeLines(message, con = file.path(outdir, "diffbind.error.txt"))
}

fail_or_record <- function(status, reason, message) {
    if (allow_partial) {
        if (identical(status, "FAIL") && !is.null(message) && nzchar(message)) {
            write_error(message)
        }
        write_summary(status, reason)
        quit(save = "no")
    }
    stop(message)
}

check_peaks <- function(paths) {
    issues <- list(missing = character(), empty = character())
    if (is.null(paths)) {
        return(issues)
    }
    for (path in paths) {
        if (is.na(path) || !nzchar(path)) {
            issues$missing <- c(issues$missing, "<missing>")
            next
        }
        if (!file.exists(path)) {
            issues$missing <- c(issues$missing, path)
            next
        }
        lines <- tryCatch(readLines(path, n = 2, warn = FALSE), error = function(e) character())
        nonempty <- lines[!grepl("^\\s*$|^#", lines)]
        if (length(nonempty) == 0) {
            issues$empty <- c(issues$empty, path)
        }
    }
    issues
}

peaks_issues <- check_peaks(samplesheet$Peaks)
if (length(peaks_issues$missing) > 0 || length(peaks_issues$empty) > 0) {
    reason <- if (length(peaks_issues$missing) > 0) "missing_peaks" else "empty_peaks"
    message_parts <- c()
    if (length(peaks_issues$missing) > 0) {
        message_parts <- c(message_parts, sprintf("missing peak files: %s", paste(unique(peaks_issues$missing), collapse = ", ")))
    }
    if (length(peaks_issues$empty) > 0) {
        message_parts <- c(message_parts, sprintf("empty peak files: %s", paste(unique(peaks_issues$empty), collapse = ", ")))
    }
    detail <- paste(message_parts, collapse = "; ")
    status <- if (allow_partial) "SKIP" else "FAIL"
    fail_or_record(status, reason, detail)
}

if (!requireNamespace("DiffBind", quietly = TRUE)) {
    fail_or_record("FAIL", "diffbind_unavailable", "DiffBind package not available.")
}

call_diffbind <- function(fun_name, args) {
    if (!exists(fun_name, envir = asNamespace("DiffBind"), inherits = FALSE)) {
        stop(sprintf("DiffBind function %s is not available.", fun_name))
    }
    fun <- get(fun_name, envir = asNamespace("DiffBind"))
    formal_names <- names(formals(fun))
    args <- args[names(args) %in% formal_names]
    do.call(fun, args)
}

plot_to_pdf <- function(filename, fun_name, args) {
    plots_dir <- file.path(outdir, "plots")
    if (!dir.exists(plots_dir)) {
        dir.create(plots_dir, recursive = TRUE)
    }
    pdf(file.path(plots_dir, filename))
    on.exit(dev.off(), add = TRUE)
    call_diffbind(fun_name, args)
}

tryCatch({
    suppressPackageStartupMessages(library(DiffBind))
    extra_params <- resolve_constants(read_extra_params(extra_params_path))

    dba_args <- utils::modifyList(list(sampleSheet = samplesheet), extra_params_for(extra_params, "dba"))
    dba_obj <- do.call(dba, dba_args)
    summit_size <- if (recenter) summits else 0
    count_args <- utils::modifyList(list(DBA = dba_obj, summits = summit_size, minOverlap = min_overlap), extra_params_for(extra_params, "dba_count"))
    dba_obj <- do.call(dba.count, count_args)

    if (use_spikein && "spikein_scale_factor" %in% colnames(records)) {
        factors <- suppressWarnings(as.numeric(records$spikein_scale_factor))
        if (any(!is.na(factors))) {
            size_factors <- ifelse(factors == 0, 1, 1 / factors)
            names(size_factors) <- records$sample_id
            size_factors <- size_factors[samplesheet$SampleID]
            norm_args <- utils::modifyList(list(DBA = dba_obj, normalize = DBA_NORM_LIB, library = size_factors), extra_params_for(extra_params, "dba_normalize"))
            dba_obj <- do.call(dba.normalize, norm_args)
            norm_out <- data.frame(sample_id = samplesheet$SampleID, size_factor = size_factors, stringsAsFactors = FALSE)
            write.table(norm_out, file = file.path(outdir, "diffbind.normalization_factors.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
        }
    } else if (!is.null(norm_method) && norm_method != "native") {
        norm_args <- utils::modifyList(list(DBA = dba_obj, method = norm_method), extra_params_for(extra_params, "dba_normalize"))
        dba_obj <- do.call(dba.normalize, norm_args)
    }

    contrast_args <- utils::modifyList(list(DBA = dba_obj, categories = DBA_CONDITION, group1 = contrast_labels[1], group2 = contrast_labels[2]), extra_params_for(extra_params, "dba_contrast"))
    dba_obj <- do.call(dba.contrast, contrast_args)

    method_flag <- ifelse(toupper(backend) == "EDGER", DBA_EDGER, DBA_DESEQ2)
    analyze_args <- utils::modifyList(list(DBA = dba_obj, method = method_flag), extra_params_for(extra_params, "dba_analyze"))
    dba_obj <- do.call(dba.analyze, analyze_args)

    plot_to_pdf("PCA.pdf", "dba.plotPCA", list(DBA = dba_obj, attributes = DBA_CONDITION, label = DBA_ID))
    plot_to_pdf("correlation_heatmap.pdf", "dba.plotHeatmap", list(DBA = dba_obj, correlations = TRUE))
    plot_to_pdf("MA.pdf", "dba.plotMA", list(DBA = dba_obj, contrast = 1))
    plot_to_pdf("volcano.pdf", "dba.plotVolcano", list(DBA = dba_obj, contrast = 1))

    report_args <- utils::modifyList(list(DBA = dba_obj, th = fdr, fold = lfc), extra_params_for(extra_params, "dba_report"))
    report <- do.call(dba.report, report_args)
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

    write_summary(
        status = "RUN",
        reason = "ok",
        n_tested = nrow(results_out),
        n_fdr_pass = nrow(sig),
        n_up = sum(sig$log2FC >= lfc),
        n_down = sum(sig$log2FC <= -lfc)
    )

    saveRDS(dba_obj, file = file.path(outdir, "diffbind.dba.rds"))

}, error = function(e) {
    fail_or_record("FAIL", "diffbind_error", e$message)
})
