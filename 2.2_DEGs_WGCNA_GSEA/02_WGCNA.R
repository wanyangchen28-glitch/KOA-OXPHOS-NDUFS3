args <- commandArgs(trailingOnly = FALSE)

script_file <- sub("^--file=", "", grep("^--file=", args, value = TRUE)[1])

code_root <- dirname(dirname(normalizePath(script_file, mustWork = TRUE)))

source(file.path(code_root, "00_configuration.R"))

suppressPackageStartupMessages(library(WGCNA))

discovery_dir <- file.path(input_root, "bulk_discovery")

outdir <- file.path(output_root, "WGCNA")

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

expr <- as.matrix(read.csv(file.path(discovery_dir, "00_preprocess", "06_ComBat_HGNC_maxIQR_expression.csv"), row.names = 1, 
    check.names = FALSE))

meta <- read.csv(file.path(discovery_dir, "group_metadata.csv"))

meta <- meta[match(colnames(expr), meta$Sample), , drop = FALSE]

stopifnot(identical(colnames(expr), meta$Sample), !anyNA(meta$Group))

datExpr <- as.data.frame(t(expr))

gene_mad <- apply(datExpr, 2, mad, na.rm = TRUE)

cutoff <- quantile(gene_mad, 0.5, na.rm = TRUE, type = 7)

datExpr <- datExpr[, is.finite(gene_mad) & gene_mad >= cutoff, drop = FALSE]

qc <- goodSamplesGenes(datExpr, verbose = 3)

datExpr <- datExpr[qc$goodSamples, qc$goodGenes, drop = FALSE]

meta <- meta[match(rownames(datExpr), meta$Sample), , drop = FALSE]

datExpr[] <- sapply(datExpr, as.numeric)

network <- blockwiseModules(datExpr, power = 12, minModuleSize = 100, deepSplit = 2, mergeCutHeight = 0.25, numericLabels = TRUE, 
    networkType = "signed", maxBlockSize = ncol(datExpr), pamRespectsDendro = FALSE, saveTOMs = FALSE, loadTOMs = FALSE, 
    verbose = 3)

colors <- labels2colors(network$colors)

MEs <- orderMEs(moduleEigengenes(datExpr, colors = colors)$eigengenes)

phenotype <- as.integer(meta$Group == "OA")

module_cor <- stats::cor(MEs, phenotype, method = "pearson")

module_p <- corPvalueStudent(module_cor, nrow(datExpr))

module_table <- data.frame(module = sub("^ME", "", rownames(module_cor)), cor = as.numeric(module_cor), P = as.numeric(module_p))

module_table$adjusted_P <- p.adjust(module_table$P, "BH")

module_table$selected <- module_table$module != "grey" & module_table$cor >= 0.7 & module_table$adjusted_P < 0.05

gene_table <- data.frame(gene = colnames(datExpr), module = colors)

membership <- stats::cor(datExpr, MEs, method = "pearson", use = "p")

significance <- as.numeric(stats::cor(datExpr, phenotype, method = "pearson", use = "p"))

gene_table$gene_significance <- significance

for (m in colnames(MEs)) gene_table[[paste0("membership_", m)]] <- membership[, m]

write.csv(module_table, file.path(outdir, "module_trait_correlations.csv"), row.names = FALSE)

write.csv(gene_table, file.path(outdir, "gene_module_membership.csv"), row.names = FALSE)

saveRDS(network, file.path(outdir, "network.rds"))

