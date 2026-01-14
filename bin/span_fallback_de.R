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
counts_path <- params[["counts"]]
samples_path <- params[["samples"]]
contrast <- params[["contrast"]]
backend <- params[["backend"]]
use_spikein <- as_bool(params[["use-spikein"]])
out_path <- params[["out"]]

if (is.null(counts_path) || is.null(samples_path) || is.null(out_path) || is.null(contrast)) {
    stop("Missing required arguments: --counts, --samples, --contrast, --out")
}

contrast_parts <- strsplit(contrast, ",")[[1]]
contrast_parts <- trimws(contrast_parts)
if (length(contrast_parts) < 2) {
    stop("Contrast must contain treated,control")
}
treated <- contrast_parts[1]
control <- contrast_parts[2]

counts <- read.delim(counts_path, check.names = FALSE, stringsAsFactors = FALSE)
if (nrow(counts) == 0) {
    write.table(data.frame(chr=character(), start=integer(), end=integer(), log2FC=numeric(), pval=numeric(), FDR=numeric()),
                file=out_path, sep="\t", quote=FALSE, row.names=FALSE)
    quit(save = "no")
}

region <- counts[, 1:3]
count_mat <- as.matrix(counts[, 4:ncol(counts), drop = FALSE])
colnames(count_mat) <- colnames(counts)[4:ncol(counts)]

samples <- read.delim(samples_path, stringsAsFactors = FALSE)
if (!"sample_id" %in% colnames(samples)) {
    stop("Samples file must contain sample_id")
}

sample_order <- match(colnames(count_mat), samples$sample_id)
if (any(is.na(sample_order))) {
    stop("Counts contain sample IDs missing from samples metadata")
}

samples <- samples[sample_order, ]
coldata <- data.frame(condition = factor(samples$condition, levels = c(control, treated)))
rownames(coldata) <- samples$sample_id

size_factors <- NULL
if (use_spikein && "spikein_scale_factor" %in% colnames(samples)) {
    factors <- suppressWarnings(as.numeric(samples$spikein_scale_factor))
    if (any(!is.na(factors))) {
        size_factors <- ifelse(factors == 0, 1, 1 / factors)
        size_factors[is.na(size_factors)] <- 1
        size_factors[!is.finite(size_factors)] <- 1
    }
}

if (is.null(backend)) {
    backend <- "DESeq2"
}

backend <- as.character(backend)

if (backend == "edgeR") {
    if (!requireNamespace("edgeR", quietly = TRUE)) {
        stop("edgeR package not available")
    }
    y <- edgeR::DGEList(counts = count_mat, group = coldata$condition)
    if (!is.null(size_factors)) {
        norm <- size_factors / exp(mean(log(size_factors)))
        y$samples$norm.factors <- norm
    } else {
        y <- edgeR::calcNormFactors(y, method = "TMM")
    }
    design <- model.matrix(~coldata$condition)
    y <- edgeR::estimateDisp(y, design)
    fit <- edgeR::glmQLFit(y, design)
    res <- edgeR::glmQLFTest(fit, coef = 2)
    tbl <- res$table
    tbl$FDR <- p.adjust(tbl$PValue, method = "BH")
    log2fc <- tbl$logFC
    pval <- tbl$PValue
    fdr <- tbl$FDR
} else {
    if (!requireNamespace("DESeq2", quietly = TRUE)) {
        stop("DESeq2 package not available")
    }
    dds <- DESeq2::DESeqDataSetFromMatrix(countData = count_mat, colData = coldata, design = ~condition)
    if (!is.null(size_factors)) {
        DESeq2::sizeFactors(dds) <- size_factors
    }
    dds <- DESeq2::DESeq(dds)
    res <- DESeq2::results(dds, contrast = c("condition", treated, control))
    log2fc <- res$log2FoldChange
    pval <- res$pvalue
    fdr <- res$padj
}

out <- data.frame(
    chr = region[[1]],
    start = region[[2]],
    end = region[[3]],
    log2FC = log2fc,
    pval = pval,
    FDR = fdr,
    stringsAsFactors = FALSE
)

write.table(out, file = out_path, sep = "\t", quote = FALSE, row.names = FALSE)
