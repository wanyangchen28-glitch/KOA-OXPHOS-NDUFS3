args <- commandArgs(trailingOnly = FALSE)

script_file <- sub("^--file=", "", grep("^--file=", args, value = TRUE)[1])

code_root <- dirname(dirname(normalizePath(script_file, mustWork = TRUE)))

source(file.path(code_root, "00_configuration.R"))

options(stringsAsFactors = FALSE, width = 160)

project_root <- input_root

dataset_dir <- file.path(project_root, "bulk_training")

old_training_dir <- file.path(project_root, "bulk_discovery")

gse32317_rma_file <- file.path(dataset_dir, "GSE32317_RAW", "GSE32317_complete_RMA_data.txt")

old_expression_file <- file.path(old_training_dir, "00_preprocess", "02_RMA_HGNC_maxIQR_expression.csv")

old_metadata_file <- file.path(old_training_dir, "00_preprocess", "04_sample_info.tsv")

expected_old_counts <- c(Control = 20L, OA = 20L)

expected_new_counts <- c(Control = 7L, OA = 19L)

expected_total_counts <- c(Control = 27L, OA = 39L)

outdir <- file.path(dataset_dir, "00_preprocess")

summary_dir <- file.path(dataset_dir, "summary")

log_dir <- file.path(summary_dir, "logs")

local_lib <- file.path(dataset_dir, "R_library")

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

.libPaths(c(local_lib, .libPaths()))

suppressPackageStartupMessages({
    library(AnnotationDbi)
    library(hgu133plus2.db)
    library(sva)
    library(ggplot2)
    library(gridExtra)
})

if (!all(file.exists(c(gse32317_rma_file, old_expression_file, old_metadata_file)))) {
    stop("Missing required input file(s). Check USER CONFIGURATION paths.")
}

read_expression <- function(path) {
    x <- read.csv(path, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
    x <- as.matrix(x)
    storage.mode(x) <- "numeric"
    x
}

read_old_metadata <- function(path) {
    x <- read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
    if (!all(c("GSM", "Group", "Batch") %in% names(x))) {
        stop("Old training metadata must contain GSM, Group, and Batch.")
    }
    x <- x[x$Group %in% c("Control", "OA"), c("GSM", "Group", "Batch"), drop = FALSE]
    names(x)[names(x) == "GSM"] <- "Sample"
    x$Platform <- "GPL96"
    x$OriginalLabel <- NA
    x
}

read_gse32317_rma <- function(path) {
    z <- read.delim(path, header = FALSE, check.names = FALSE, quote = "", comment.char = "", stringsAsFactors = FALSE)
    if (nrow(z) < 3L || ncol(z) != 27L) {
        stop("Unexpected GSE32317 supplementary-matrix dimensions.")
    }
    conditions <- as.character(unlist(z[1, -1, drop = TRUE]))
    samples <- as.character(unlist(z[2, -1, drop = TRUE]))
    probes <- as.character(z[-c(1, 2), 1])
    expr <- as.matrix(z[-c(1, 2), -1, drop = FALSE])
    storage.mode(expr) <- "numeric"
    rownames(expr) <- probes
    colnames(expr) <- samples
    group <- ifelse(grepl("^OA_", conditions), "OA", ifelse(grepl("^Healthy_", conditions), "Control", NA))
    meta <- data.frame(Sample = samples, Group = group, Batch = "GSE32317", Platform = "GPL570", OriginalLabel = conditions, 
        stringsAsFactors = FALSE)
    if (anyNA(meta$Group) || anyDuplicated(meta$Sample)) {
        stop("GSE32317 group parsing or sample identifiers failed.")
    }
    list(expr_probe = expr, metadata = meta)
}

collapse_gpl570_to_symbol <- function(expr_probe) {
    probe_ids <- rownames(expr_probe)
    ann <- AnnotationDbi::select(hgu133plus2.db, keys = probe_ids, columns = "SYMBOL", keytype = "PROBEID")
    ann <- ann[!is.na(ann$SYMBOL) & nzchar(ann$SYMBOL), c("PROBEID", "SYMBOL"), drop = FALSE]
    ann <- ann[!duplicated(ann), , drop = FALSE]
    n_symbol <- table(ann$PROBEID)
    ann <- ann[ann$PROBEID %in% names(n_symbol[n_symbol == 1L]), , drop = FALSE]
    ann <- ann[match(intersect(probe_ids, ann$PROBEID), ann$PROBEID), , drop = FALSE]
    mapping_rate <- nrow(ann)/length(probe_ids)
    if (mapping_rate < 0.3) 
        stop("GPL570 probe-to-HGNC mapping rate < 30%.")
    ann$IQR <- apply(expr_probe[ann$PROBEID, , drop = FALSE], 1, IQR, na.rm = TRUE)
    ann <- ann[order(ann$SYMBOL, -ann$IQR, ann$PROBEID), , drop = FALSE]
    ann_max <- ann[!duplicated(ann$SYMBOL), , drop = FALSE]
    expr_gene <- expr_probe[ann_max$PROBEID, , drop = FALSE]
    rownames(expr_gene) <- ann_max$SYMBOL
    list(expression = expr_gene, annotation = ann_max, mapping_rate = mapping_rate)
}

save_plot <- function(plot, stem, width, height) {
    ggsave(paste0(stem, ".pdf"), plot, width = width, height = height, units = "in")
    ggsave(paste0(stem, ".png"), plot, width = width, height = height, units = "in", dpi = 300)
}

pca_data <- function(x, metadata, stage) {
    pc <- prcomp(t(x), center = TRUE, scale. = FALSE)
    explained <- 100 * pc$sdev^2/sum(pc$sdev^2)
    out <- data.frame(Sample = rownames(pc$x), PC1 = pc$x[, 1], PC2 = pc$x[, 2], Group = metadata$Group[match(rownames(pc$x), 
        metadata$Sample)], Batch = metadata$Batch[match(rownames(pc$x), metadata$Sample)], Stage = stage, stringsAsFactors = FALSE)
    attr(out, "explained") <- explained
    out
}

plot_pca <- function(z, title) {
    explained <- attr(z, "explained")
    ggplot(z, aes(PC1, PC2, colour = Group, shape = Batch)) + geom_point(size = 2.8, alpha = 0.9) + scale_colour_manual(values = c(Control = "#2C7FB8", 
        OA = "#D95F02")) + labs(title = title, x = sprintf("PC1 (%.1f%%)", explained[1]), y = sprintf("PC2 (%.1f%%)", explained[2])) + 
        theme_classic(base_size = 11) + theme(plot.title = element_text(hjust = 0.5, face = "bold"))
}

factor_p <- function(z, variable) {
    c(PC1_p = summary(aov(z$PC1 ~ z[[variable]]))[[1]][1, "Pr(>F)"], PC2_p = summary(aov(z$PC2 ~ z[[variable]]))[[1]][1, 
        "Pr(>F)"])
}

old_expr <- read_expression(old_expression_file)

old_meta <- read_old_metadata(old_metadata_file)

old_meta <- old_meta[match(colnames(old_expr), old_meta$Sample), , drop = FALSE]

if (!identical(old_meta$Sample, colnames(old_expr))) stop("Old expression/metadata sample mismatch.")

if (!identical(as.integer(table(old_meta$Group)[names(expected_old_counts)]), as.integer(expected_old_counts))) {
    stop("Unexpected old training group counts.")
}

gse32317 <- read_gse32317_rma(gse32317_rma_file)

new_meta <- gse32317$metadata

if (!identical(as.integer(table(new_meta$Group)[names(expected_new_counts)]), as.integer(expected_new_counts))) {
    stop("Unexpected GSE32317 group counts.")
}

new_processed <- collapse_gpl570_to_symbol(gse32317$expr_probe)

new_expr <- new_processed$expression

new_meta <- new_meta[match(colnames(new_expr), new_meta$Sample), , drop = FALSE]

if (!identical(new_meta$Sample, colnames(new_expr))) stop("GSE32317 expression/metadata sample mismatch.")

write.csv(gse32317$expr_probe, file.path(outdir, "01_GSE32317_supplied_RMA_probe_expression.csv"), quote = FALSE)

write.csv(new_expr, file.path(outdir, "02_GSE32317_HGNC_maxIQR_expression.csv"), quote = FALSE)

write.table(new_processed$annotation, file.path(outdir, "03_GSE32317_probe_to_HGNC_maxIQR.tsv"), sep = "\t", quote = FALSE, 
    row.names = FALSE)

write.table(new_meta, file.path(outdir, "04_GSE32317_sample_info.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

common_genes <- intersect(rownames(old_expr), rownames(new_expr))

if (length(common_genes) < 10000L) stop("Too few common HGNC genes across GPL96/GPL570: ", length(common_genes))

expr_precombat <- cbind(old_expr[common_genes, old_meta$Sample, drop = FALSE], new_expr[common_genes, new_meta$Sample, drop = FALSE])

metadata <- rbind(old_meta, new_meta)

metadata <- metadata[match(colnames(expr_precombat), metadata$Sample), , drop = FALSE]

if (!identical(metadata$Sample, colnames(expr_precombat))) stop("Combined expression/metadata mismatch.")

if (anyDuplicated(metadata$Sample)) stop("Duplicated sample IDs in expanded training cohort.")

if (!identical(as.integer(table(metadata$Group)[names(expected_total_counts)]), as.integer(expected_total_counts))) {
    stop("Unexpected total group counts after expansion.")
}

batch <- factor(metadata$Batch, levels = c("GSE55235", "GSE55457", "GSE32317"))

group <- factor(metadata$Group, levels = c("Control", "OA"))

mod <- model.matrix(~group)

expr_combat <- sva::ComBat(dat = expr_precombat, batch = batch, mod = mod, par.prior = TRUE, prior.plots = FALSE, mean.only = FALSE)

colnames(expr_combat) <- colnames(expr_precombat)

rownames(expr_combat) <- rownames(expr_precombat)

write.csv(expr_precombat, file.path(outdir, "05_expanded_training_preComBat_HGNC_expression.csv"), quote = FALSE)

write.csv(expr_combat, file.path(outdir, "06_expanded_training_ComBat_HGNC_expression.csv"), quote = FALSE)

write.csv(expr_combat, file.path(outdir, "expression_matrix.csv"), quote = FALSE)

write.table(metadata, file.path(outdir, "07_expanded_training_sample_info.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

write.csv(metadata[, c("Sample", "Group", "Batch", "Platform")], file.path(outdir, "group_metadata.csv"), row.names = FALSE)

write.table(metadata, file.path(dataset_dir, "sample_info.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

module_counts <- as.data.frame.matrix(table(metadata$Batch, metadata$Group))

module_counts$Batch <- rownames(module_counts)

module_counts <- module_counts[, c("Batch", "Control", "OA")]

write.table(module_counts, file.path(summary_dir, "batch_group_counts.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

pca_before <- pca_data(expr_precombat, metadata, "Before ComBat")

pca_after <- pca_data(expr_combat, metadata, "After ComBat")

pca_all <- rbind(pca_before, pca_after)

write.table(pca_all, file.path(outdir, "08_PCA_scores_before_after.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

pca_qc <- data.frame(Stage = c("Before ComBat", "After ComBat"), PC1_batch_p = c(factor_p(pca_before, "Batch")["PC1_p"], 
    factor_p(pca_after, "Batch")["PC1_p"]), PC2_batch_p = c(factor_p(pca_before, "Batch")["PC2_p"], factor_p(pca_after, "Batch")["PC2_p"]), 
    PC1_group_p = c(factor_p(pca_before, "Group")["PC1_p"], factor_p(pca_after, "Group")["PC1_p"]), PC2_group_p = c(factor_p(pca_before, 
        "Group")["PC2_p"], factor_p(pca_after, "Group")["PC2_p"]), row.names = NULL)

write.table(pca_qc, file.path(outdir, "09_PCA_batch_group_QC.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

p_before <- plot_pca(pca_before, "Expanded training cohort before ComBat")

p_after <- plot_pca(pca_after, "Expanded training cohort after ComBat")

p_pair <- gridExtra::arrangeGrob(p_before, p_after, ncol = 2)

ggsave(file.path(outdir, "10_PCA_before_after_ComBat.pdf"), p_pair, width = 11, height = 5.2, units = "in")

ggsave(file.path(outdir, "10_PCA_before_after_ComBat.png"), p_pair, width = 11, height = 5.2, units = "in", dpi = 300)

distribution_summary <- function(x, stage) {
    do.call(rbind, lapply(seq_len(ncol(x)), function(i) {
        data.frame(Stage = stage, Sample = colnames(x)[i], Group = metadata$Group[match(colnames(x)[i], metadata$Sample)], 
            Batch = metadata$Batch[match(colnames(x)[i], metadata$Sample)], Median = median(x[, i]), IQR = IQR(x[, i]), stringsAsFactors = FALSE)
    }))
}

distribution <- rbind(distribution_summary(expr_precombat, "Before ComBat"), distribution_summary(expr_combat, "After ComBat"))

write.table(distribution, file.path(outdir, "11_expression_distribution_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

invisible(NULL)

invisible(NULL)

invisible(NULL)

invisible(NULL)

capture.output(sessionInfo(), file = file.path(summary_dir, "sessionInfo.txt"))

message("Completed GSE32317-expanded training preprocessing: ", outdir)

