#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) %% 2 != 0) {
  stop("Arguments must be provided as --key value pairs")
}
arg_map <- list()
for (i in seq(1, length(args), by = 2)) {
  key <- sub("^--", "", args[i])
  arg_map[[key]] <- args[i + 1]
}

matrix_path <- arg_map[["matrix"]]
samplesheet_path <- arg_map[["samplesheet"]]
treated_label <- arg_map[["treated"]]
control_label <- arg_map[["control"]]
bootstrap <- as.integer(arg_map[["bootstrap"]])
k_value <- as.integer(arg_map[["k-value"]])
out_path <- arg_map[["out"]]

if (is.null(matrix_path) || is.null(samplesheet_path) || is.null(out_path)) {
  stop("Missing required arguments")
}

suppressMessages(library(ROTS))

mat <- read.table(matrix_path, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
mat <- as.matrix(mat)

samples <- read.csv(samplesheet_path, header = TRUE, stringsAsFactors = FALSE)
subset <- samples[samples$condition %in% c(treated_label, control_label), ]
if (nrow(subset) == 0) {
  stop("No samples match contrast labels")
}

sample_ids <- subset$sample_id
missing <- setdiff(sample_ids, colnames(mat))
if (length(missing) > 0) {
  stop(paste("Missing samples in matrix:", paste(missing, collapse = ",")))
}

mat <- mat[, sample_ids, drop = FALSE]

classes <- ifelse(subset$condition == treated_label, 1, 2)

rot <- ROTS(mat, groups = classes, B = bootstrap, K = k_value)

pvals <- NULL
fdrs <- NULL

if (!is.null(rot$pvalue)) {
  pvals <- rot$pvalue
} else if (!is.null(rot$pval)) {
  pvals <- rot$pval
}

if (is.null(pvals)) {
  if ("pvalue" %in% slotNames(rot)) {
    pvals <- slot(rot, "pvalue")
  }
}

if (is.null(pvals)) {
  if (exists("pvalue", where = asNamespace("ROTS"), mode = "function")) {
    pvals <- ROTS::pvalue(rot)
  }
}

if (!is.null(rot$FDR)) {
  fdrs <- rot$FDR
} else if (!is.null(rot$fdr)) {
  fdrs <- rot$fdr
}

if (is.null(fdrs) && !is.null(pvals)) {
  fdrs <- p.adjust(pvals, method = "BH")
}

if (is.null(pvals)) {
  stop("Unable to extract p-values from ROTS output")
}

out <- data.frame(
  bin_id = rownames(mat),
  pval = as.numeric(pvals),
  FDR = as.numeric(fdrs)
)

write.table(out, file = out_path, sep = "\t", quote = FALSE, row.names = FALSE)
