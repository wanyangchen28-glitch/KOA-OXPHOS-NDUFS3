args <- commandArgs(trailingOnly = FALSE)

script_file <- sub("^--file=", "", grep("^--file=", args, value = TRUE)[1])

code_root <- dirname(dirname(normalizePath(script_file, mustWork = TRUE)))

source(file.path(code_root, "00_configuration.R"))

invisible(gc())

suppressPackageStartupMessages({
    library(Seurat)
    library(ggplot2)
})

root_dir <- file.path(input_root, "singlecell")

object_file <- file.path(root_dir, "08_PC30_res0.6_rebuild_DotPlot/01_rds", "GSE216651_scRNA_PC30_res0.6_6major_before_UMAP_adjustment.rds")

umap_file <- file.path(root_dir, "10_PC30_res0.6_adjusted_UMAP/01_rds", "PC30_original_UMAP_embeddings_seed20260716.rds")

out_dir <- file.path(root_dir, "12_remove_cluster16_17_and_replot")

rds_dir <- file.path(out_dir, "01_rds")

table_dir <- file.path(out_dir, "02_tables")

fig_dir <- file.path(out_dir, "03_figures")

log_dir <- file.path(out_dir, "04_logs")

for (d in c(rds_dir, table_dir, fig_dir, log_dir)) {
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

if (!file.exists(object_file)) stop("未找到Step 08对象：\n", object_file)

if (!file.exists(umap_file)) stop("未找到原始UMAP坐标：\n", umap_file)

message("读取Step 08对象，请耐心等待……")

obj <- readRDS(object_file)

if (!inherits(obj, "Seurat")) stop("输入文件不是Seurat对象。")

cluster_candidates <- intersect(c("PC30_snn_res.0.6", "cluster_res_0_6", "integrated_snn_res.0.6", "seurat_clusters"), colnames(obj@meta.data))

if (length(cluster_candidates) == 0) {
    stop("找不到resolution 0.6对应的cluster列。")
}

cluster_col <- cluster_candidates[1]

major_candidates <- intersect(c("major_cell_type", "major_celltype", "cell_type_6", "celltype_6", "cell_type", "celltype"), 
    colnames(obj@meta.data))

if (length(major_candidates) == 0) stop("找不到六大细胞类型注释列。")

major_col <- major_candidates[1]

cluster_vec <- as.character(obj@meta.data[[cluster_col]])

remove_clusters <- c("16", "17")

if (!all(remove_clusters %in% unique(cluster_vec))) {
    stop("对象中没有同时找到cluster 16和17。实际cluster：\n", paste(sort(unique(cluster_vec)), collapse = ", "))
}

remove_cells <- colnames(obj)[cluster_vec %in% remove_clusters]

keep_cells <- setdiff(colnames(obj), remove_cells)

removed_metadata_columns <- intersect(c(cluster_col, major_col, "sample_id", "donor_id", "Group", "nCount_RNA", "nFeature_RNA", 
    "percent.mt", "scDblFinder.score", "scDblFinder.class"), colnames(obj@meta.data))

removed_metadata <- obj@meta.data[remove_cells, removed_metadata_columns, drop = FALSE]

removed_metadata$cell_barcode <- rownames(removed_metadata)

removed_metadata <- removed_metadata[, c("cell_barcode", removed_metadata_columns), drop = FALSE]

write.csv(removed_metadata, file.path(table_dir, "removed_cluster16_17_cell_metadata.csv"), row.names = FALSE)

writeLines(keep_cells, file.path(table_dir, "retained_cell_barcodes.txt"))

removal_summary <- data.frame(item = c("Cells before removal", "Cluster 16 removed", "Cluster 17 removed", "Total removed", 
    "Cells retained"), n_cells = c(ncol(obj), sum(cluster_vec == "16"), sum(cluster_vec == "17"), length(remove_cells), length(keep_cells)), 
    stringsAsFactors = FALSE)

removal_summary$fraction_of_original <- removal_summary$n_cells/ncol(obj)

write.csv(removal_summary, file.path(table_dir, "cluster16_17_removal_summary.csv"), row.names = FALSE)

print(removal_summary, row.names = FALSE)

obj_clean <- subset(obj, cells = keep_cells)

obj_clean@meta.data[[cluster_col]] <- droplevels(factor(as.character(obj_clean@meta.data[[cluster_col]])))

umap_in <- readRDS(umap_file)

if (inherits(umap_in, "DimReduc")) {
    umap_mat <- Embeddings(umap_in)
} else {
    umap_mat <- as.matrix(umap_in)
}

if (is.null(rownames(umap_mat))) stop("UMAP坐标没有细胞barcode行名。")

if (!all(colnames(obj_clean) %in% rownames(umap_mat))) {
    stop("UMAP坐标与过滤后对象的细胞barcode不一致。")
}

umap_mat <- umap_mat[colnames(obj_clean), 1:2, drop = FALSE]

colnames(umap_mat) <- c("UMAP_1", "UMAP_2")

obj_clean[["umap_original_filtered"]] <- CreateDimReducObject(embeddings = umap_mat, key = "UMAPF_", assay = DefaultAssay(obj_clean))

major_levels <- c("Fibroblasts", "Endothelial", "Mural", "Myeloid", "Lymphoid", "Mast")

obj_clean@meta.data[[major_col]] <- factor(as.character(obj_clean@meta.data[[major_col]]), levels = major_levels)

major_colors <- c(Fibroblasts = "#4E79A7", Endothelial = "#59A14F", Mural = "#F28E2B", Myeloid = "#E15759", Lymphoid = "#B07AA1", 
    Mast = "#76B7B2")

set.seed(20260716)

p_umap <- DimPlot(obj_clean, reduction = "umap_original_filtered", group.by = major_col, cols = major_colors, label = FALSE, 
    raster = TRUE, raster.dpi = c(900, 900), pt.size = 0.1) + labs(x = "UMAP 1", y = "UMAP 2", colour = "Cell type") + coord_equal() + 
    theme_classic(base_size = 8, base_family = "sans") + theme(axis.line = element_line(linewidth = 0.45, colour = "black"), 
    axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_text(size = 9, colour = "black"), legend.title = element_text(size = 8.5), 
    legend.text = element_text(size = 8), legend.key.height = grid::unit(3.5, "mm"), plot.margin = margin(4, 4, 4, 4, unit = "mm"))

final_pdf <- file.path(fig_dir, "UMAP_6_major_celltypes_PC30_after_removing_clusters16_17.pdf")

ggsave(filename = final_pdf, plot = p_umap, width = 135, height = 105, units = "mm", device = grDevices::cairo_pdf, bg = "white")

filtered_rds <- file.path(rds_dir, "GSE216651_scRNA_PC30_res0.6_without_clusters16_17.rds")

rm(obj)

invisible(gc())

message("正在保存过滤后的Seurat对象；该文件较大，请耐心等待……")

saveRDS(obj_clean, filtered_rds, compress = "gzip")

stopifnot(ncol(obj_clean) == 60302, !any(as.character(obj_clean@meta.data[[cluster_col]]) %in% remove_clusters), file.exists(final_pdf), 
    file.exists(filtered_rds))

invisible(NULL)

sink(file.path(log_dir, "sessionInfo.txt"))

print(sessionInfo())

sink()

message("\nStep 12完成。")

message("过滤后对象：", filtered_rds)

message("新UMAP图：", final_pdf)

