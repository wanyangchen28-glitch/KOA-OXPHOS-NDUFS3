args <- commandArgs(trailingOnly = FALSE)

script_file <- sub("^--file=", "", grep("^--file=", args, value = TRUE)[1])

code_root <- dirname(dirname(normalizePath(script_file, mustWork = TRUE)))

source(file.path(code_root, "00_configuration.R"))

suppressPackageStartupMessages({
    library(Seurat)
    library(dplyr)
})

atlas <- readRDS(file.path(input_root, "singlecell", "GSE216651_scRNA_PC30_res0.6_without_clusters16_17.rds"))

score_file <- file.path(output_root, "33_phosphorylation_score_major_and_fibroblast", "01_whole_atlas_phosphorylation_score_5methods.csv")

scores <- read.csv(score_file)

genes <- c("ATP6V1A", "NDUFS3")

expr <- LayerData(atlas, assay = "SCT", layer = "data")[genes, , drop = FALSE]

meta <- atlas[[]]

meta$cell_id <- rownames(meta)

dat <- inner_join(meta, scores[, c("cell_id", "score_consensus_z")], by = "cell_id")

dat$score_group <- ifelse(dat$score_consensus_z > 0, "score_UP", "score_DOWN")

for (g in genes) dat[[g]] <- as.numeric(expr[g, dat$cell_id])

outdir <- file.path(output_root, "candidate_gene_expression")

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

write.csv(dat[, c("cell_id", "Group", "major_cell_type", "score_consensus_z", "score_group", genes)], file.path(outdir, "all_cell_gene_expression.csv"), 
    row.names = FALSE)

fib <- readRDS(file.path(input_root, "singlecell", "GSE216651_Fibroblasts_final_7_published_subtypes.rds"))

fib_expr <- LayerData(fib, assay = "SCT", layer = "data")[genes, , drop = FALSE]

fib_meta <- fib[[]]

fib_meta$cell_id <- rownames(fib_meta)

for (g in genes) fib_meta[[g]] <- as.numeric(fib_expr[g, fib_meta$cell_id])

labels <- read.csv(file.path(output_root, "Scissor", "cell_labels.csv"))

fib_meta <- left_join(fib_meta, labels[, c("cell_id", "Scissor_raw_class")], by = "cell_id")

prg4 <- fib_meta[fib_meta$fibroblast_subtype_published == "PRG4+ lining fibroblasts" & fib_meta$Scissor_raw_class %in% c("Scissor+", 
    "Scissor-"), ]

tests <- do.call(rbind, lapply(genes, function(g) {
    p <- wilcox.test(prg4[[g]][prg4$Scissor_raw_class == "Scissor+"], prg4[[g]][prg4$Scissor_raw_class == "Scissor-"], exact = FALSE, 
        alternative = "two.sided")$p.value
    data.frame(Gene = g, P = p)
}))

tests$adjusted_P <- p.adjust(tests$P, "BH")

write.csv(fib_meta[, c("cell_id", "Group", "fibroblast_subtype_published", "Scissor_raw_class", genes)], file.path(outdir, 
    "fibroblast_gene_expression.csv"), row.names = FALSE)

write.csv(tests, file.path(outdir, "PRG4_Scissor_expression_statistics.csv"), row.names = FALSE)

