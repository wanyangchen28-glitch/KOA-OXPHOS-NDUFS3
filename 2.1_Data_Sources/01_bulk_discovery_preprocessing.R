args <- commandArgs(trailingOnly = FALSE)

script_file <- sub("^--file=", "", grep("^--file=", args, value = TRUE)[1])

code_root <- dirname(dirname(normalizePath(script_file, mustWork = TRUE)))

source(file.path(code_root, "00_configuration.R"))

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
    library(affy)
    library(Biobase)
    library(AnnotationDbi)
    library(hgu133a.db)
    library(sva)
    library(ggplot2)
    library(gridExtra)
})

root <- file.path(input_root, "bulk_discovery")

out <- file.path(root, "00_preprocess")

dir.create(out, recursive = TRUE, showWarnings = FALSE)

read_pheno <- function(path, gse, group_col) {
    x <- read.delim(path, row.names = 1, check.names = FALSE, quote = "", comment.char = "")
    raw <- tolower(as.character(x[[group_col]]))
    group <- ifelse(grepl("osteoarthrit", raw), "OA", ifelse(grepl("healthy|normal.*control", raw), "Control", NA))
    keep <- !is.na(group)
    if (!any(keep)) 
        stop("No OA/control samples found in ", gse)
    data.frame(GSM = rownames(x)[keep], Group = group[keep], Batch = gse, OriginalLabel = as.character(x[[group_col]][keep]), 
        stringsAsFactors = FALSE, row.names = rownames(x)[keep])
}

p1 <- read_pheno(file.path(root, "GSE55235_RAW", "GSE55235_pheno.tsv"), "GSE55235", "disease state:ch1")

p2 <- read_pheno(file.path(root, "GSE55457_RAW", "GSE55457_pheno.tsv"), "GSE55457", "clinical status:ch1")

sample_info <- rbind(p1, p2)

if (anyDuplicated(sample_info$GSM)) stop("Duplicated GSM IDs in metadata")

if (nrow(sample_info) != 40L || any(table(sample_info$Group) == 0)) stop("Expected 40 OA/control samples with nonzero groups")

cel_files <- unlist(lapply(c("GSE55235", "GSE55457"), function(g) {
    list.files(file.path(root, paste0(g, "_RAW")), pattern = "(?i)\\.cel$", full.names = TRUE)
}))

cel_gsm <- sub("_.*$", "", basename(cel_files))

cel_files <- cel_files[cel_gsm %in% sample_info$GSM]

cel_gsm <- sub("_.*$", "", basename(cel_files))

if (length(cel_files) != 40L || anyDuplicated(cel_gsm)) stop("CEL/metadata mismatch")

cel_files <- cel_files[match(sample_info$GSM, cel_gsm)]

message("Reading ", length(cel_files), " CEL files and applying joint RMA...")

raw <- ReadAffy(filenames = cel_files)

sampleNames(raw) <- sample_info$GSM

rma_eset <- rma(raw, background = TRUE, normalize = TRUE, pmcorrect.method = "pmonly", normalize.type = "quantiles", summary.method = "medianpolish")

expr_probe <- exprs(rma_eset)

colnames(expr_probe) <- sample_info$GSM

rownames(expr_probe) <- sub("\\.at$", "", rownames(expr_probe))

probe_ids <- rownames(expr_probe)

ann <- AnnotationDbi::select(hgu133a.db, keys = probe_ids, columns = "SYMBOL", keytype = "PROBEID")

ann <- ann[!is.na(ann$SYMBOL) & nzchar(ann$SYMBOL), c("PROBEID", "SYMBOL")]

ann <- ann[!duplicated(ann), ]

sym_n <- table(ann$PROBEID)

ann <- ann[ann$PROBEID %in% names(sym_n[sym_n == 1L]), ]

ann <- ann[match(intersect(probe_ids, ann$PROBEID), ann$PROBEID), ]

if (nrow(ann)/length(probe_ids) < 0.3) stop("HGNC/SYMBOL mapping rate < 30%")

expr_annot <- expr_probe[ann$PROBEID, , drop = FALSE]

expr_by_symbol <- expr_annot

rownames(expr_by_symbol) <- ann$SYMBOL

ann$IQR <- apply(expr_by_symbol, 1, IQR, na.rm = TRUE)

ann <- ann[order(ann$SYMBOL, -ann$IQR, ann$PROBEID), ]

ann_max <- ann[!duplicated(ann$SYMBOL), ]

expr_gene <- expr_probe[ann_max$PROBEID, , drop = FALSE]

rownames(expr_gene) <- ann_max$SYMBOL

sample_info <- sample_info[match(colnames(expr_gene), sample_info$GSM), ]

rownames(sample_info) <- sample_info$GSM

stopifnot(identical(rownames(sample_info), colnames(expr_gene)))

write.csv(expr_probe, file.path(out, "01_RMA_probe_expression.csv"), quote = FALSE)

write.csv(expr_gene, file.path(out, "02_RMA_HGNC_maxIQR_expression.csv"), quote = FALSE)

write.table(ann_max, file.path(out, "03_probe_to_HGNC_maxIQR.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

write.table(sample_info, file.path(out, "04_sample_info.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

saveRDS(rma_eset, file.path(out, "05_joint_RMA_ExpressionSet.rds"))

batch <- factor(sample_info$Batch)

group <- factor(sample_info$Group, levels = c("Control", "OA"))

expr_combat <- sva::ComBat(dat = as.matrix(expr_gene), batch = batch, mod = model.matrix(~group), par.prior = TRUE, prior.plots = FALSE, 
    mean.only = FALSE)

colnames(expr_combat) <- colnames(expr_gene)

rownames(expr_combat) <- rownames(expr_gene)

write.csv(expr_combat, file.path(out, "06_ComBat_HGNC_maxIQR_expression.csv"), quote = FALSE)

pca_scores <- function(x, stage) {
    pc <- prcomp(t(x), center = TRUE, scale. = FALSE)
    ve <- 100 * pc$sdev^2/sum(pc$sdev^2)
    z <- data.frame(Sample = rownames(pc$x), PC1 = pc$x[, 1], PC2 = pc$x[, 2], Group = sample_info[rownames(pc$x), "Group"], 
        Batch = sample_info[rownames(pc$x), "Batch"], Stage = stage)
    attr(z, "variance") <- ve
    z
}

plot_pca <- function(z) {
    ve <- attr(z, "variance")
    ggplot(z, aes(PC1, PC2, color = Group, shape = Batch)) + geom_point(size = 3.2) + labs(x = sprintf("PC1 (%.1f%%)", ve[1]), 
        y = sprintf("PC2 (%.1f%%)", ve[2])) + theme_classic(base_size = 12)
}

pc_before <- pca_scores(expr_gene, "Before ComBat")

pc_after <- pca_scores(expr_combat, "After ComBat")

pc <- rbind(pc_before, pc_after)

write.table(pc, file.path(out, "07_PCA_scores_before_after.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

pdf(file.path(out, "08_PCA_before_after_ComBat.pdf"), width = 11, height = 5.5)

grid.arrange(plot_pca(pc_before), plot_pca(pc_after), ncol = 2)

dev.off()

ggsave(file.path(out, "08_PCA_before_after_ComBat.png"), arrangeGrob(plot_pca(pc_before), plot_pca(pc_after), ncol = 2), 
    width = 11, height = 5.5, dpi = 300)

factor_p <- function(z, factor_name) {
    f <- z[[factor_name]]
    c(PC1_p = summary(aov(PC1 ~ f, data = z))[[1]]["f", "Pr(>F)"], PC2_p = summary(aov(PC2 ~ f, data = z))[[1]]["f", "Pr(>F)"])
}

qc <- data.frame(Stage = c("Before ComBat", "After ComBat"), PC1_batch_p = c(factor_p(pc_before, "Batch")["PC1_p"], factor_p(pc_after, 
    "Batch")["PC1_p"]), PC2_batch_p = c(factor_p(pc_before, "Batch")["PC2_p"], factor_p(pc_after, "Batch")["PC2_p"]), PC1_group_p = c(factor_p(pc_before, 
    "Group")["PC1_p"], factor_p(pc_after, "Group")["PC1_p"]), PC2_group_p = c(factor_p(pc_before, "Group")["PC2_p"], factor_p(pc_after, 
    "Group")["PC2_p"]), check.names = FALSE)

write.table(qc, file.path(out, "09_PCA_batch_separation_QC.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

make_distribution_long <- function(x, stage) {
    do.call(rbind, lapply(seq_len(ncol(x)), function(i) {
        data.frame(Stage = stage, Sample = colnames(x)[i], Group = sample_info[colnames(x)[i], "Group"], Batch = sample_info[colnames(x)[i], 
            "Batch"], Expression = as.numeric(x[, i]), stringsAsFactors = FALSE)
    }))
}

dist_long <- rbind(make_distribution_long(expr_gene, "Before ComBat"), make_distribution_long(expr_combat, "After ComBat"))

dist_long$Stage <- factor(dist_long$Stage, levels = c("Before ComBat", "After ComBat"))

dist_summary <- do.call(rbind, lapply(split(dist_long, dist_long$Stage), function(d) {
    do.call(rbind, lapply(split(d, d$Sample), function(s) {
        data.frame(Stage = as.character(s$Stage[1]), Sample = s$Sample[1], Group = s$Group[1], Batch = s$Batch[1], Median = median(s$Expression), 
            IQR = IQR(s$Expression))
    }))
}))

write.table(dist_summary, file.path(out, "12_expression_distribution_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

p_box <- ggplot(dist_long, aes(x = Sample, y = Expression, fill = Group)) + geom_boxplot(outlier.size = 0.08, linewidth = 0.15, 
    width = 0.78) + facet_grid(. ~ Stage, scales = "free_x", space = "free_x") + scale_fill_manual(values = c(Control = "#F8766D", 
    OA = "#00BFC4")) + labs(x = NULL, y = "RMA expression") + theme_classic(base_size = 9) + theme(axis.text.x = element_blank(), 
    axis.ticks.x = element_blank(), legend.position = "top", strip.background = element_blank(), strip.text = element_text(face = "bold"))

make_density <- function(x, stage) {
    do.call(rbind, lapply(seq_len(ncol(x)), function(i) {
        d <- density(as.numeric(x[, i]), na.rm = TRUE, n = 256)
        data.frame(Stage = stage, Sample = colnames(x)[i], x = d$x, y = d$y, Group = sample_info[colnames(x)[i], "Group"], 
            Batch = sample_info[colnames(x)[i], "Batch"])
    }))
}

density_long <- rbind(make_density(expr_gene, "Before ComBat"), make_density(expr_combat, "After ComBat"))

density_long$Stage <- factor(density_long$Stage, levels = c("Before ComBat", "After ComBat"))

p_density <- ggplot(density_long, aes(x, y, group = Sample, color = Batch)) + geom_line(alpha = 0.3, linewidth = 0.25) + 
    facet_grid(. ~ Stage) + scale_color_manual(values = c(GSE55235 = "#303030", GSE55457 = "#2C7FB8")) + labs(x = "RMA expression", 
    y = "Density", color = "Batch") + theme_classic(base_size = 9) + theme(legend.position = "top", strip.background = element_blank(), 
    strip.text = element_text(face = "bold"))

pdf(file.path(out, "13_expression_distribution_before_after_ComBat.pdf"), width = 11, height = 8)

grid.arrange(p_box, p_density, ncol = 1, heights = c(1.15, 1))

dev.off()

ggsave(file.path(out, "13_expression_distribution_before_after_ComBat.png"), arrangeGrob(p_box, p_density, ncol = 1, heights = c(1.15, 
    1)), width = 11, height = 8, dpi = 300)

writeLines(capture.output(sessionInfo()), file.path(out, "10_sessionInfo.txt"))

