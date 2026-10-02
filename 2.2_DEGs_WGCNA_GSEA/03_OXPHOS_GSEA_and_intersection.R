args <- commandArgs(trailingOnly = FALSE)

script_file <- sub("^--file=", "", grep("^--file=", args, value = TRUE)[1])

code_root <- dirname(dirname(normalizePath(script_file, mustWork = TRUE)))

source(file.path(code_root, "00_configuration.R"))

suppressPackageStartupMessages({
    library(limma)
    library(fgsea)
})

discovery_dir <- file.path(input_root, "bulk_discovery")

outdir <- file.path(output_root, "OXPHOS_GSEA")

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

expr <- as.matrix(read.csv(file.path(discovery_dir, "00_preprocess", "06_ComBat_HGNC_maxIQR_expression.csv"), row.names = 1, 
    check.names = FALSE))

meta <- read.delim(file.path(discovery_dir, "00_preprocess", "04_sample_info.tsv"))

meta <- meta[match(colnames(expr), meta$GSM), , drop = FALSE]

stopifnot(identical(colnames(expr), meta$GSM), !anyNA(meta$Group))

group <- factor(meta$Group, levels = c("Control", "OA"))

design <- model.matrix(~0 + group)

colnames(design) <- levels(group)

fit <- eBayes(contrasts.fit(lmFit(expr, design), makeContrasts(OA - Control, levels = design)))

tt <- topTable(fit, coef = 1, number = Inf, sort.by = "none")

tt$gene <- toupper(trimws(rownames(tt)))

ranks <- setNames(tt$t, tt$gene)

ranks <- sort(ranks[is.finite(ranks)], decreasing = TRUE)

gmt <- file.path(discovery_dir, "08_Hub_gene_GSEA", "00_gene_sets", "KEGG_official_canonical_hsa_20260717.symbols.gmt")

lines <- strsplit(readLines(gmt), "\t", fixed = TRUE)

line <- Filter(function(x) x[1] == "KEGG_00190_OXIDATIVE_PHOSPHORYLATION", lines)

stopifnot(length(line) == 1)

genes <- unique(toupper(trimws(line[[1]][-(1:2)])))

set.seed(20260802)

gsea <- fgseaMultilevel(list(OXPHOS = intersect(genes, rownames(expr))), stats = ranks, minSize = 10, maxSize = 1000, eps = 0)

modules <- read.csv(file.path(output_root, "WGCNA", "module_trait_correlations.csv"))

module_genes <- read.csv(file.path(output_root, "WGCNA", "gene_module_membership.csv"))

selected_modules <- modules[modules$selected, , drop = FALSE]
selected_module <- selected_modules$module[which.max(selected_modules$cor)]
selected_genes <- module_genes$gene[module_genes$module == selected_module]

up <- tt$gene[tt$logFC >= 0.5 & tt$adj.P.Val < 0.05]

selected <- Reduce(intersect, list(genes, up, selected_genes))

candidates <- tt[match(selected, tt$gene), c("gene", "logFC", "adj.P.Val")]

names(candidates) <- c("gene", "bulk_log2FC", "bulk_FDR")

candidates$leading_edge <- candidates$gene %in% gsea$leadingEdge[[1]]

candidates$WGCNA_modules <- module_genes$module[match(candidates$gene, module_genes$gene)]

gsea_table <- as.data.frame(gsea)

gsea_table$leadingEdge <- vapply(gsea_table$leadingEdge, paste, character(1), collapse = ";")

write.csv(gsea_table, file.path(outdir, "OXPHOS_GSEA.csv"), row.names = FALSE)

write.csv(candidates, file.path(outdir, "candidate_genes.csv"), row.names = FALSE)

write.csv(tt, file.path(outdir, "limma_gene_statistics.csv"), row.names = FALSE)

