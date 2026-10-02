args <- commandArgs(trailingOnly = FALSE)

script_file <- sub("^--file=", "", grep("^--file=", args, value = TRUE)[1])

code_root <- dirname(dirname(normalizePath(script_file, mustWork = TRUE)))

source(file.path(code_root, "00_configuration.R"))

invisible(gc())

ROOT_DIR <- file.path(input_root, "singlecell")

INPUT_RDS <- file.path(ROOT_DIR, "20_Fibroblast_secondary_clustering", "01_resolution_screen", "01_rds", "GSE216651_Fibroblasts_PC20_resolution_0.2_to_0.6_screen.rds")

OUT_DIR <- file.path(ROOT_DIR, "20_Fibroblast_secondary_clustering", "03_clean_recluster_after_cluster8_removal_corrected")

RDS_DIR <- file.path(OUT_DIR, "01_rds")

TABLE_DIR <- file.path(OUT_DIR, "02_tables")

FIGURE_DIR <- file.path(OUT_DIR, "03_figures")

LOG_DIR <- file.path(OUT_DIR, "04_logs")

OLD_CLUSTER_FIELD <- "fibro_res_0_5"

REMOVE_CLUSTER <- "8"

PRE_CLEAN_CLUSTER_FIELD <- "preclean_fibro_res_0_5"

FINAL_CLUSTER_FIELD <- "fibro_clean_res_0_5"

CLUSTER_ASSAY <- "integrated"

PCA_REDUCTION <- "fibro_clean_pca50"

UMAP_REDUCTION <- "fibro_clean_umap_PC20"

PCA_DIMS <- 1:20

N_PCS_TO_CALCULATE <- 50L

GRAPH_NN <- "fibro_clean_PC20_nn"

GRAPH_SNN <- "fibro_clean_PC20_snn"

K_NEIGHBORS <- 20L

FINAL_RESOLUTION <- 0.5

CLUSTER_ALGORITHM <- 1L

UMAP_N_NEIGHBORS <- 30L

UMAP_MIN_DIST <- 0.3

UMAP_METRIC <- "cosine"

SAMPLE_FIELD <- "sample_id"

GROUP_FIELD <- "Group"

RANDOM_SEED <- 20260730L

EXPECTED_INPUT_CELLS <- 32531L

EXPECTED_REMOVED_CELLS <- 1527L

EXPECTED_FINAL_CELLS <- 31004L

CLUSTER_PALETTE <- c("#D55E00", "#0072B2", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#8A2BE2", "#7F3C8D", "#11A579", "#3969AC", 
    "#F2B701", "#E73F74", "#80BA5A", "#E68310", "#008695")

SAMPLE_PALETTE <- c("#D55E00", "#0072B2", "#009E73", "#CC79A7", "#E69F00", "#56B4E9")

required_packages <- c("Seurat", "SeuratObject", "ggplot2", "future")

missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing_packages) > 0L) {
    stop("缺少 R 包：", paste(missing_packages, collapse = ", "), "。请先安装后重新运行。")
}

suppressPackageStartupMessages({
    library(Seurat)
    library(ggplot2)
})

if (utils::packageVersion("Seurat") < "5.0.0") {
    stop("本脚本按 Seurat v5 编写；当前版本为 ", packageVersion("Seurat"), "。")
}

for (d in c(RDS_DIR, TABLE_DIR, FIGURE_DIR, LOG_DIR)) {
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

numeric_cluster_levels <- function(x) {
    x <- unique(as.character(x))
    if (all(grepl("^[0-9]+$", x))) {
        as.character(sort(as.integer(x)))
    }
    else {
        sort(x)
    }
}

save_pdf_png <- function(plot, stem, width_mm = 145, height_mm = 115, dpi = 600) {
    ggsave(filename = file.path(FIGURE_DIR, paste0(stem, ".pdf")), plot = plot, width = width_mm, height = height_mm, units = "mm", 
        device = grDevices::cairo_pdf, bg = "white")
    ggsave(filename = file.path(FIGURE_DIR, paste0(stem, ".png")), plot = plot, width = width_mm, height = height_mm, units = "mm", 
        dpi = dpi, bg = "white")
}

theme_umap <- function() {
    theme_classic(base_size = 10, base_family = "sans") + theme(axis.line = element_line(linewidth = 0.5, colour = "black"), 
        axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_text(size = 10, colour = "black"), 
        legend.title = element_text(size = 9), legend.text = element_text(size = 8.5), plot.title = element_text(size = 12, 
            face = "bold", hjust = 0.5), plot.subtitle = element_text(size = 9, colour = "#3F3F3F", hjust = 0.5), plot.margin = margin(5, 
            5, 5, 5, unit = "mm"))
}

run_analysis <- function() {
    log_file <- file.path(LOG_DIR, "20E_remove_cluster8_recluster_run.log")
    log_con <- file(log_file, open = "wt")
    sink(log_con, type = "output", split = TRUE)
    sink(log_con, type = "message")
    on.exit({
        while (sink.number(type = "message") > 2L) sink(type = "message")
        while (sink.number(type = "output") > 0L) sink(type = "output")
        close(log_con)
    }, add = TRUE)
    cat("开始时间：", format(Sys.time()), "\n", sep = "")
    cat("Seurat：", as.character(packageVersion("Seurat")), "\n", sep = "")
    cat("输入对象：", INPUT_RDS, "\n", sep = "")
    cat("删除条件：", OLD_CLUSTER_FIELD, " == ", REMOVE_CLUSTER, "\n", sep = "")
    cat("固定参数：dims 1:20；resolution 0.5\n\n")
    if (!file.exists(INPUT_RDS)) 
        stop("未找到输入 RDS：\n", INPUT_RDS)
    message("读取 resolution screen 对象……")
    fibro <- readRDS(INPUT_RDS)
    if (!inherits(fibro, "Seurat")) 
        stop("输入文件不是 Seurat 对象。")
    if (ncol(fibro) != EXPECTED_INPUT_CELLS) {
        stop("输入对象为 ", ncol(fibro), " cells，与预期的 ", EXPECTED_INPUT_CELLS, " 不一致。")
    }
    required_meta <- c(OLD_CLUSTER_FIELD, SAMPLE_FIELD, GROUP_FIELD)
    missing_meta <- setdiff(required_meta, colnames(fibro@meta.data))
    if (length(missing_meta) > 0L) {
        stop("对象缺少必要 metadata：", paste(missing_meta, collapse = ", "))
    }
    if (!CLUSTER_ASSAY %in% SeuratObject::Assays(fibro)) {
        stop("对象中不存在 assay：", CLUSTER_ASSAY)
    }
    if (!"scale.data" %in% SeuratObject::Layers(fibro[[CLUSTER_ASSAY]])) {
        stop("integrated assay 缺少既有 scale.data，不能按20A一致方案重算 PCA。")
    }
    old_clusters <- as.character(fibro@meta.data[[OLD_CLUSTER_FIELD]])
    names(old_clusters) <- colnames(fibro)
    remove_cells <- names(old_clusters)[old_clusters == REMOVE_CLUSTER]
    keep_cells <- names(old_clusters)[old_clusters != REMOVE_CLUSTER]
    if (length(remove_cells) != EXPECTED_REMOVED_CELLS) {
        stop("识别到的 cluster 8 为 ", length(remove_cells), " cells，与预期的 ", EXPECTED_REMOVED_CELLS, " 不一致；为防止误删，分析停止。")
    }
    if (length(keep_cells) != EXPECTED_FINAL_CELLS) {
        stop("删除后应为 ", EXPECTED_FINAL_CELLS, " cells，实际为 ", length(keep_cells), "。")
    }
    fibro@meta.data[[PRE_CLEAN_CLUSTER_FIELD]] <- old_clusters
    message("删除 cluster 8：", length(remove_cells), " cells。")
    fibro_clean <- subset(fibro, cells = keep_cells)
    removed_meta <- fibro@meta.data[remove_cells, , drop = FALSE]
    input_sample_counts <- table(as.character(fibro@meta.data[[SAMPLE_FIELD]]))
    removed_sample_counts <- table(as.character(removed_meta[[SAMPLE_FIELD]]))
    sample_levels <- sort(unique(as.character(fibro@meta.data[[SAMPLE_FIELD]])))
    removed_by_sample <- data.frame(sample_id = sample_levels, input_fibro_cells = as.integer(input_sample_counts[sample_levels]), 
        removed_cluster8_cells = as.integer(removed_sample_counts[sample_levels]), stringsAsFactors = FALSE)
    removed_by_sample$removed_fraction_within_sample <- removed_by_sample$removed_cluster8_cells/removed_by_sample$input_fibro_cells
    sample_group_map <- unique(fibro@meta.data[, c(SAMPLE_FIELD, GROUP_FIELD), drop = FALSE])
    if (anyDuplicated(sample_group_map[[SAMPLE_FIELD]]) > 0L) {
        stop("同一个 sample_id 对应多个 Group，metadata 不一致。")
    }
    removed_by_sample <- merge(removed_by_sample, sample_group_map, by.x = "sample_id", by.y = SAMPLE_FIELD, all.x = TRUE, 
        sort = FALSE)
    write.csv(removed_by_sample, file.path(TABLE_DIR, "01_removed_cluster8_by_sample.csv"), row.names = FALSE)
    write.csv(data.frame(cell_barcode = remove_cells, sample_id = as.character(removed_meta[[SAMPLE_FIELD]]), Group = as.character(removed_meta[[GROUP_FIELD]]), 
        removed_cluster = REMOVE_CLUSTER, stringsAsFactors = FALSE), file.path(TABLE_DIR, "02_removed_cluster8_cell_barcodes.csv"), 
        row.names = FALSE)
    rm(fibro, removed_meta)
    invisible(gc())
    old_resolution_fields <- grep("^fibro_res_[0-9]+_[0-9]+$", colnames(fibro_clean@meta.data), value = TRUE)
    for (field in old_resolution_fields) fibro_clean[[field]] <- NULL
    old_reductions <- SeuratObject::Reductions(fibro_clean)
    for (reduction_name in old_reductions) {
        fibro_clean[[reduction_name]] <- NULL
    }
    fibro_clean@graphs <- list()
    fibro_clean@neighbors <- list()
    DefaultAssay(fibro_clean) <- CLUSTER_ASSAY
    pca_features <- SeuratObject::VariableFeatures(fibro_clean[[CLUSTER_ASSAY]])
    if (length(pca_features) == 0L) {
        pca_features <- rownames(fibro_clean[[CLUSTER_ASSAY]])
    }
    pca_features <- intersect(pca_features, rownames(fibro_clean[[CLUSTER_ASSAY]]))
    if (length(pca_features) < 500L) {
        stop("可用于 PCA 的 integrated features 少于500个，分析停止。")
    }
    if (length(pca_features) <= N_PCS_TO_CALCULATE) {
        stop("PCA features 数量不足以稳定计算50 PCs。")
    }
    future::plan(future::sequential)
    options(future.globals.maxSize = 16 * 1024^3)
    set.seed(RANDOM_SEED)
    message("沿用既有 integrated scale.data，重新计算 PCA50（不重新 ScaleData）……")
    fibro_clean <- RunPCA(object = fibro_clean, assay = CLUSTER_ASSAY, features = pca_features, npcs = N_PCS_TO_CALCULATE, 
        reduction.name = PCA_REDUCTION, reduction.key = "fcPC_", seed.use = RANDOM_SEED, verbose = TRUE)
    if (ncol(Embeddings(fibro_clean, PCA_REDUCTION)) < max(PCA_DIMS)) {
        stop("重新计算的 PCA 少于20个维度。")
    }
    message("使用固定 PC1:20 重建邻接图……")
    fibro_clean <- FindNeighbors(object = fibro_clean, reduction = PCA_REDUCTION, dims = PCA_DIMS, k.param = K_NEIGHBORS, 
        graph.name = c(GRAPH_NN, GRAPH_SNN), verbose = TRUE)
    set.seed(RANDOM_SEED)
    message("使用固定 PC1:20 重算 UMAP……")
    fibro_clean <- RunUMAP(object = fibro_clean, reduction = PCA_REDUCTION, dims = PCA_DIMS, n.neighbors = UMAP_N_NEIGHBORS, 
        min.dist = UMAP_MIN_DIST, metric = UMAP_METRIC, umap.method = "uwot", reduction.name = UMAP_REDUCTION, reduction.key = "fcUMAP20_", 
        seed.use = RANDOM_SEED, verbose = TRUE)
    set.seed(RANDOM_SEED)
    message("运行一次聚类：resolution = 0.5……")
    fibro_clean <- FindClusters(object = fibro_clean, graph.name = GRAPH_SNN, resolution = FINAL_RESOLUTION, algorithm = CLUSTER_ALGORITHM, 
        random.seed = RANDOM_SEED, cluster.name = FINAL_CLUSTER_FIELD, verbose = TRUE)
    final_levels <- numeric_cluster_levels(fibro_clean@meta.data[[FINAL_CLUSTER_FIELD]])
    fibro_clean@meta.data[[FINAL_CLUSTER_FIELD]] <- factor(as.character(fibro_clean@meta.data[[FINAL_CLUSTER_FIELD]]), levels = final_levels)
    Idents(fibro_clean) <- FINAL_CLUSTER_FIELD
    if (ncol(fibro_clean) != EXPECTED_FINAL_CELLS) {
        stop("重聚类对象细胞数发生异常变化。")
    }
    final_cluster <- as.character(fibro_clean@meta.data[[FINAL_CLUSTER_FIELD]])
    final_counts <- as.data.frame(table(cluster = factor(final_cluster, levels = final_levels)), stringsAsFactors = FALSE)
    names(final_counts)[2] <- "n_cells"
    final_counts$fraction_of_clean_fibroblasts <- final_counts$n_cells/ncol(fibro_clean)
    write.csv(final_counts, file.path(TABLE_DIR, "03_final_cluster_cell_counts.csv"), row.names = FALSE)
    by_sample <- as.data.frame(table(cluster = factor(final_cluster, levels = final_levels), sample_id = as.character(fibro_clean@meta.data[[SAMPLE_FIELD]])), 
        stringsAsFactors = FALSE)
    names(by_sample)[3] <- "n_cells"
    cluster_totals <- setNames(final_counts$n_cells, final_counts$cluster)
    by_sample$within_cluster_fraction <- by_sample$n_cells/cluster_totals[as.character(by_sample$cluster)]
    by_sample <- merge(by_sample, sample_group_map, by.x = "sample_id", by.y = SAMPLE_FIELD, all.x = TRUE, sort = FALSE)
    write.csv(by_sample, file.path(TABLE_DIR, "04_final_cluster_by_sample_counts.csv"), row.names = FALSE)
    transition <- as.data.frame(table(preclean_cluster = as.character(fibro_clean@meta.data[[PRE_CLEAN_CLUSTER_FIELD]]), 
        final_cluster = final_cluster), stringsAsFactors = FALSE)
    names(transition)[3] <- "n_cells"
    transition <- transition[transition$n_cells > 0L, , drop = FALSE]
    preclean_totals <- tapply(transition$n_cells, transition$preclean_cluster, sum)
    transition$fraction_of_preclean_cluster <- transition$n_cells/preclean_totals[transition$preclean_cluster]
    write.csv(transition, file.path(TABLE_DIR, "05_preclean_to_final_cluster_transition.csv"), row.names = FALSE)
    cluster_colours <- setNames(CLUSTER_PALETTE[seq_along(final_levels)], final_levels)
    p_cluster <- DimPlot(object = fibro_clean, reduction = UMAP_REDUCTION, group.by = FINAL_CLUSTER_FIELD, cols = cluster_colours, 
        label = TRUE, repel = TRUE, label.size = 4.3, label.box = TRUE, label.color = "black", raster = TRUE, raster.dpi = c(300, 
            300), pt.size = 0.45, alpha = 0.95) + labs(title = "Fibroblasts after removal of low-quality cluster 8", subtitle = "PC1:20; resolution 0.5; n = 31,004", 
        x = "UMAP 1", y = "UMAP 2", colour = "Cluster") + coord_equal() + theme_umap()
    save_pdf_png(p_cluster, "01_clean_UMAP_final_clusters")
    final_sample_levels <- sort(unique(as.character(fibro_clean@meta.data[[SAMPLE_FIELD]])))
    if (length(final_sample_levels) > length(SAMPLE_PALETTE)) {
        stop("样本数超过预设颜色数量。")
    }
    sample_colours <- setNames(SAMPLE_PALETTE[seq_along(final_sample_levels)], final_sample_levels)
    p_sample <- DimPlot(object = fibro_clean, reduction = UMAP_REDUCTION, group.by = SAMPLE_FIELD, cols = sample_colours, 
        shuffle = TRUE, seed = RANDOM_SEED, raster = TRUE, raster.dpi = c(300, 300), pt.size = 0.42, alpha = 0.92) + labs(title = "Clean Fibroblast UMAP by sample", 
        subtitle = "Technical distribution check", x = "UMAP 1", y = "UMAP 2", colour = "Sample") + coord_equal() + theme_umap()
    save_pdf_png(p_sample, "02_clean_UMAP_by_sample")
    group_levels <- unique(as.character(fibro_clean@meta.data[[GROUP_FIELD]]))
    if (!setequal(group_levels, c("Control", "OA"))) {
        stop("Group 应仅包含 Control 和 OA；实际为：", paste(sort(group_levels), collapse = ", "))
    }
    p_group <- DimPlot(object = fibro_clean, reduction = UMAP_REDUCTION, group.by = GROUP_FIELD, cols = c(Control = "#4E79A7", 
        OA = "#E15759"), shuffle = TRUE, seed = RANDOM_SEED, raster = TRUE, raster.dpi = c(300, 300), pt.size = 0.42, alpha = 0.92) + 
        labs(title = "Clean Fibroblast UMAP by group", subtitle = "Descriptive visualization only", x = "UMAP 1", y = "UMAP 2", 
            colour = "Group") + coord_equal() + theme_umap()
    save_pdf_png(p_group, "03_clean_UMAP_by_Group")
    output_rds <- file.path(RDS_DIR, "GSE216651_Fibroblasts_without_cluster8_PC20_res0.5.rds")
    message("保存最终清洁重聚类对象……")
    saveRDS(fibro_clean, output_rds, compress = FALSE)
    invisible(NULL)
    invisible(NULL)
    writeLines(capture.output(sessionInfo()), file.path(LOG_DIR, "20E_sessionInfo.txt"), useBytes = TRUE)
    cat("\n完成时间：", format(Sys.time()), "\n", sep = "")
    cat("删除细胞：", length(remove_cells), "\n", sep = "")
    cat("最终细胞：", ncol(fibro_clean), "\n", sep = "")
    cat("最终 cluster 数：", length(final_levels), "\n", sep = "")
    cat("最终 RDS：", output_rds, "\n", sep = "")
    cat("结果目录：", OUT_DIR, "\n", sep = "")
    cat("本脚本未进行正式亚群注释或候选基因分析。\n")
    invisible(fibro_clean)
}

result <- run_analysis()

