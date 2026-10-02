args <- commandArgs(trailingOnly = FALSE)

script_file <- sub("^--file=", "", grep("^--file=", args, value = TRUE)[1])

code_root <- dirname(dirname(normalizePath(script_file, mustWork = TRUE)))

source(file.path(code_root, "00_configuration.R"))

gc()

invisible(NULL)

library(Seurat)

library(SeuratObject)

library(Matrix)

library(dplyr)

library(tidyr)

library(ggplot2)

library(AUCell)

library(UCell)

library(singscore)

library(GSVA)

library(clusterProfiler)

atlas_file <- file.path(input_root, "singlecell/GSE216651_scRNA_PC30_res0.6_without_clusters16_17.rds")

fibro_file <- file.path(input_root, "singlecell/GSE216651_Fibroblasts_final_7_published_subtypes.rds")

phosphorylation_gmt <- file.path(input_root, "singlecell/h.all.v2026.1.Hs.symbols.gmt.txt")

outdir <- file.path(output_root, "33_phosphorylation_score_major_and_fibroblast")

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

assay_name <- "SCT"

layer_name <- "data"

group_field <- "Group"

major_celltype_field <- "major_cell_type"

fibro_subtype_field <- "fibroblast_subtype_published"

phosphorylation <- clusterProfiler::read.gmt(phosphorylation_gmt)

phosphorylation_genes <- phosphorylation %>% dplyr::filter(term == "HALLMARK_OXIDATIVE_PHOSPHORYLATION") %>% dplyr::pull(gene) %>% 
    unique()

score_five_methods <- function(object, object_name) {
    DefaultAssay(object) <- assay_name
    expr_mat <- LayerData(object = object, assay = assay_name, layer = layer_name)
    sig_genes <- intersect(phosphorylation_genes, rownames(expr_mat))
    gene_sets <- list(PHOSPHORYLATION = sig_genes)
    cat(object_name, " retained genes:", length(sig_genes), "\n")
    cells_rankings <- AUCell_buildRankings(expr_mat, nCores = 1, plotStats = TRUE)
    cells_AUC <- AUCell_calcAUC(gene_sets, cells_rankings, nCores = 1, aucMaxRank = nrow(expr_mat) * 0.1)
    object$score_AUCell <- as.numeric(getAUC(cells_AUC)["PHOSPHORYLATION", ])
    rm(cells_rankings, cells_AUC)
    gc()
    object <- AddModuleScore_UCell(object, features = gene_sets, name = NULL)
    object$score_UCell <- object$PHOSPHORYLATION
    object$PHOSPHORYLATION <- NULL
    expr_mat_dense <- as.matrix(expr_mat)
    ranked_expr <- rankGenes(expr_mat_dense)
    ss <- simpleScore(ranked_expr, upSet = sig_genes)
    if ("TotalScore" %in% colnames(ss)) {
        object$score_singscore <- ss$TotalScore
    }
    else {
        object$score_singscore <- ss[, 1]
    }
    rm(expr_mat_dense, ranked_expr, ss)
    gc()
    expr_mat_dense2 <- as.matrix(expr_mat)
    ssgsea_par <- ssgseaParam(exprData = expr_mat_dense2, geneSets = gene_sets)
    ssgsea_es <- gsva(ssgsea_par)
    object$score_ssGSEA <- as.numeric(ssgsea_es["PHOSPHORYLATION", ])
    rm(expr_mat_dense2, ssgsea_par, ssgsea_es)
    gc()
    object <- AddModuleScore(object = object, features = list(sig_genes), assay = assay_name, name = "phosphorylation_AddModuleScore")
    object$score_AddModuleScore <- object$phosphorylation_AddModuleScore1
    object$phosphorylation_AddModuleScore1 <- NULL
    zscore01 <- function(x) {
        as.numeric(scale(x))
    }
    minmax01 <- function(x) {
        rng <- range(x, na.rm = TRUE)
        if (diff(rng) == 0) 
            return(rep(0, length(x)))
        (x - rng[1])/diff(rng)
    }
    method_names <- c("score_AUCell", "score_UCell", "score_singscore", "score_ssGSEA", "score_AddModuleScore")
    for (method in method_names) {
        object@meta.data[[paste0(method, "_z")]] <- zscore01(object@meta.data[[method]])
        object@meta.data[[paste0(method, "_mm")]] <- minmax01(object@meta.data[[method]])
    }
    object$score_consensus_z <- rowMeans(object@meta.data[, paste0(method_names, "_z"), drop = FALSE], na.rm = TRUE)
    object$score_consensus_mm <- rowMeans(object@meta.data[, paste0(method_names, "_mm"), drop = FALSE], na.rm = TRUE)
    object
}

cat("Reading whole atlas...\n")

scRNA_atlas <- readRDS(atlas_file)

scRNA_atlas <- score_five_methods(scRNA_atlas, "Whole atlas")

cat("Reading final fibroblast object...\n")

scRNA_fibro <- readRDS(fibro_file)

scRNA_fibro <- score_five_methods(scRNA_fibro, "Final fibroblasts")

method_names <- c("score_AUCell", "score_UCell", "score_singscore", "score_ssGSEA", "score_AddModuleScore")

make_score_table <- function(object, group_field, celltype_field) {
    object@meta.data %>% tibble::rownames_to_column("cell_id") %>% dplyr::select(cell_id, all_of(group_field), all_of(celltype_field), 
        all_of(method_names), all_of(paste0(method_names, "_z")), all_of(paste0(method_names, "_mm")), score_consensus_z, 
        score_consensus_mm)
}

atlas_score_df <- make_score_table(scRNA_atlas, group_field, major_celltype_field)

fibro_score_df <- make_score_table(scRNA_fibro, group_field, fibro_subtype_field)

write.csv(atlas_score_df, file.path(outdir, "01_whole_atlas_phosphorylation_score_5methods.csv"), row.names = FALSE)

write.csv(fibro_score_df, file.path(outdir, "02_final_fibroblast_phosphorylation_score_5methods.csv"), row.names = FALSE)

score_long_atlas <- atlas_score_df %>% dplyr::select(all_of(group_field), all_of(major_celltype_field), score_consensus_z)

score_long_fibro <- fibro_score_df %>% dplyr::select(all_of(group_field), all_of(fibro_subtype_field), score_consensus_z)

atlas_summary <- score_long_atlas %>% group_by(.data[[major_celltype_field]]) %>% summarise(n_cells = n(), mean_score = mean(score_consensus_z, 
    na.rm = TRUE), median_score = median(score_consensus_z, na.rm = TRUE), Q1_score = quantile(score_consensus_z, 0.25, na.rm = TRUE), 
    Q3_score = quantile(score_consensus_z, 0.75, na.rm = TRUE), .groups = "drop")

atlas_group_summary <- score_long_atlas %>% group_by(.data[[group_field]], .data[[major_celltype_field]]) %>% summarise(n_cells = n(), 
    mean_score = mean(score_consensus_z, na.rm = TRUE), median_score = median(score_consensus_z, na.rm = TRUE), Q1_score = quantile(score_consensus_z, 
        0.25, na.rm = TRUE), Q3_score = quantile(score_consensus_z, 0.75, na.rm = TRUE), .groups = "drop")

fibro_summary <- score_long_fibro %>% group_by(.data[[fibro_subtype_field]]) %>% summarise(n_cells = n(), mean_score = mean(score_consensus_z, 
    na.rm = TRUE), median_score = median(score_consensus_z, na.rm = TRUE), Q1_score = quantile(score_consensus_z, 0.25, na.rm = TRUE), 
    Q3_score = quantile(score_consensus_z, 0.75, na.rm = TRUE), .groups = "drop")

fibro_group_summary <- score_long_fibro %>% group_by(.data[[group_field]], .data[[fibro_subtype_field]]) %>% summarise(n_cells = n(), 
    mean_score = mean(score_consensus_z, na.rm = TRUE), median_score = median(score_consensus_z, na.rm = TRUE), Q1_score = quantile(score_consensus_z, 
        0.25, na.rm = TRUE), Q3_score = quantile(score_consensus_z, 0.75, na.rm = TRUE), .groups = "drop")

write.csv(atlas_summary, file.path(outdir, "03_whole_atlas_score_summary_by_major_cell_type.csv"), row.names = FALSE)

write.csv(atlas_group_summary, file.path(outdir, "04_whole_atlas_score_summary_by_group_and_major_cell_type.csv"), row.names = FALSE)

write.csv(fibro_summary, file.path(outdir, "05_final_fibroblast_score_summary_by_subtype.csv"), row.names = FALSE)

write.csv(fibro_group_summary, file.path(outdir, "06_final_fibroblast_score_summary_by_group_and_subtype.csv"), row.names = FALSE)

atlas_umap <- FeaturePlot(scRNA_atlas, features = "score_consensus_z", reduction = "umap_original_filtered", cols = c("#F7F7F7", 
    "#FDB863", "#B2182B"), min.cutoff = "q05", max.cutoff = "q95") + ggtitle("Phosphorylation consensus score in whole atlas") + 
    theme_classic(base_size = 14)

ggsave(file.path(outdir, "07_whole_atlas_phosphorylation_consensus_UMAP.pdf"), atlas_umap, width = 8, height = 6)

ggsave(file.path(outdir, "07_whole_atlas_phosphorylation_consensus_UMAP.png"), atlas_umap, width = 8, height = 6, dpi = 300)

atlas_celltype_plot <- ggplot(atlas_score_df, aes(x = .data[[major_celltype_field]], y = score_consensus_z, fill = .data[[major_celltype_field]])) + 
    geom_boxplot(alpha = 0.85, outlier.shape = NA, width = 0.65) + theme_classic(base_size = 13) + theme(axis.text.x = element_text(angle = 45, 
    hjust = 1), legend.position = "none") + labs(x = NULL, y = "Phosphorylation consensus score (z)")

ggsave(file.path(outdir, "08_whole_atlas_phosphorylation_score_by_major_cell_type.pdf"), atlas_celltype_plot, width = 9, 
    height = 6)

ggsave(file.path(outdir, "08_whole_atlas_phosphorylation_score_by_major_cell_type.png"), atlas_celltype_plot, width = 9, 
    height = 6, dpi = 300)

fibro_umap <- FeaturePlot(scRNA_fibro, features = "score_consensus_z", reduction = "fibro_clean_umap_PC20", cols = c("#F7F7F7", 
    "#FDB863", "#B2182B"), min.cutoff = "q05", max.cutoff = "q95") + ggtitle("Phosphorylation consensus score in final fibroblasts") + 
    theme_classic(base_size = 14)

ggsave(file.path(outdir, "09_final_fibroblast_phosphorylation_consensus_UMAP.pdf"), fibro_umap, width = 8, height = 6)

ggsave(file.path(outdir, "09_final_fibroblast_phosphorylation_consensus_UMAP.png"), fibro_umap, width = 8, height = 6, dpi = 300)

fibro_subtype_plot <- ggplot(fibro_score_df, aes(x = .data[[fibro_subtype_field]], y = score_consensus_z, fill = .data[[fibro_subtype_field]])) + 
    geom_boxplot(alpha = 0.85, outlier.shape = NA, width = 0.65) + theme_classic(base_size = 13) + theme(axis.text.x = element_text(angle = 45, 
    hjust = 1), legend.position = "none") + labs(x = NULL, y = "Phosphorylation consensus score (z)")

ggsave(file.path(outdir, "10_final_fibroblast_phosphorylation_score_by_subtype.pdf"), fibro_subtype_plot, width = 10, height = 6)

ggsave(file.path(outdir, "10_final_fibroblast_phosphorylation_score_by_subtype.png"), fibro_subtype_plot, width = 10, height = 6, 
    dpi = 300)

saveRDS(scRNA_atlas, file.path(outdir, "whole_atlas_phosphorylation_scored.rds"))

saveRDS(scRNA_fibro, file.path(outdir, "final_fibroblast_phosphorylation_scored.rds"))

write.csv(data.frame(gene_set = "PHOSPHORYLATION", genes_in_gmt = length(phosphorylation_genes), genes_in_whole_atlas = length(intersect(phosphorylation_genes, 
    rownames(LayerData(scRNA_atlas[[assay_name]], layer = layer_name)))), genes_in_final_fibroblasts = length(intersect(phosphorylation_genes, 
    rownames(LayerData(scRNA_fibro[[assay_name]], layer = layer_name)))), scoring_methods = paste(method_names, collapse = ";"), 
    stringsAsFactors = FALSE), file.path(outdir, "11_gene_set_and_method_information.csv"), row.names = FALSE)

cat("Phosphorylation scoring completed.\n")

cat("Output directory:", outdir, "\n")

