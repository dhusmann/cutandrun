#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(optparse)
    library(LOLA)
    library(GenomicRanges)
    library(rtracklayer)
    library(ggplot2)
})

option_list <- list(
    make_option(c("--bed"), type = "character"),
    make_option(c("--universe"), type = "character"),
    make_option(c("--db"), type = "character"),
    make_option(c("--label"), type = "character"),
    make_option(c("--out_tsv"), type = "character"),
    make_option(c("--out_pdf"), type = "character", default = "")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (!dir.exists(opt$db)) {
    stop(paste("LOLA database not found:", opt$db))
}

bed_lines <- readLines(opt$bed, warn = FALSE)
bed_lines <- bed_lines[!grepl("^#", bed_lines)]
if (length(bed_lines) == 0) {
    write.table(data.frame(), opt$out_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
    if (opt$out_pdf != "") {
        pdf(opt$out_pdf)
        plot.new()
        text(0.5, 0.5, "No regions available for enrichment")
        dev.off()
    }
    quit(status = 0)
}

user_set <- rtracklayer::import(opt$bed, format = "bed")
universe <- rtracklayer::import(opt$universe, format = "bed")
region_db <- LOLA::loadRegionDB(opt$db)

user_sets <- GRangesList()
user_sets[[opt$label]] <- user_set

res <- LOLA::runLOLA(user_sets, universe, region_db)
if (!is.null(res) && nrow(res) > 0) {
    res <- res[res$userSet == opt$label, , drop = FALSE]
}

write.table(res, opt$out_tsv, sep = "\t", quote = FALSE, row.names = FALSE)

if (opt$out_pdf != "" && !is.null(res) && nrow(res) > 0) {
    top <- res[order(res$pValue), , drop = FALSE]
    if (nrow(top) > 20) {
        top <- top[1:20, , drop = FALSE]
    }
    top$neglog10p <- -log10(top$pValue)
    top$label <- if ("description" %in% colnames(top)) top$description else top$filename
    pdf(opt$out_pdf)
    ggplot(top, aes(x = reorder(label, neglog10p), y = neglog10p)) +
        geom_col(fill = "#4C72B0") +
        coord_flip() +
        xlab("Region set") +
        ylab("-log10(p-value)") +
        ggtitle(paste("LOLA enrichment:", opt$label))
    dev.off()
}
