args <- commandArgs(trailingOnly = FALSE)

script_file <- sub("^--file=", "", grep("^--file=", args, value = TRUE)[1])

code_root <- dirname(dirname(normalizePath(script_file, mustWork = TRUE)))

source(file.path(code_root, "00_configuration.R"))

gc()

work_dir <- file.path(output_root, "43_NDUFS3_OA_PRG4_lining_virtual_knockout")

input_rds <- file.path(input_root, "singlecell/GSE216651_Fibroblasts_final_7_published_subtypes.rds")

whole_atlas_rds <- file.path(input_root, "singlecell/GSE216651_scRNA_PC30_res0.6_without_clusters16_17.rds")

dir.create(work_dir, showWarnings = FALSE, recursive = TRUE)

setwd(work_dir)

out_dirs <- c("01_object", "02_count_matrix", "03_scTenifoldKnk_result", "04_diffRegulation_table", "05_KO_network_plot", 
    "06_volcano_plot", "07_top_gene_plot", "08_manifold_plot", "09_marker_response_plot", "10_GO_enrichment", "11_summary_tables")

for (d in out_dirs) {
    dir.create(d, showWarnings = FALSE, recursive = TRUE)
}

required_pkgs <- c("Seurat", "tidyverse", "Matrix", "ggplot2", "ggrepel", "patchwork", "pheatmap", "scTenifoldKnk", "clusterProfiler", 
    "org.Hs.eg.db", "enrichplot", "AnnotationDbi", "igraph")

missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing_pkgs) > 0) {
    stop("以下 R 包尚未安装，请先安装后再运行：\n", paste(missing_pkgs, collapse = ", "))
}

suppressPackageStartupMessages({
    library(Seurat)
    library(tidyverse)
    library(Matrix)
    library(ggplot2)
    library(ggrepel)
    library(patchwork)
    library(pheatmap)
    library(scTenifoldKnk)
    library(clusterProfiler)
    library(org.Hs.eg.db)
    library(enrichplot)
    library(AnnotationDbi)
    library(igraph)
})

assign("select", dplyr::select, envir = .GlobalEnv)

assign("filter", dplyr::filter, envir = .GlobalEnv)

assign("mutate", dplyr::mutate, envir = .GlobalEnv)

assign("arrange", dplyr::arrange, envir = .GlobalEnv)

assign("summarise", dplyr::summarise, envir = .GlobalEnv)

assign("summarize", dplyr::summarise, envir = .GlobalEnv)

assign("group_by", dplyr::group_by, envir = .GlobalEnv)

assign("left_join", dplyr::left_join, envir = .GlobalEnv)

assign("rename", dplyr::rename, envir = .GlobalEnv)

assign("E", igraph::E, envir = .GlobalEnv)

assign("V", igraph::V, envir = .GlobalEnv)

target_gene <- "NDUFS3"

locked_fibroblast_subtypes <- c("PRG4+ lining fibroblasts", "APOD+ fibroblasts", "CD34+ sublining fibroblasts", "CXCL12+ sublining fibroblasts", 
    "DKK3+ sublining fibroblasts", "RSPO3+ fibroblasts", "POSTN+ fibroblasts")

target_subtypes <- "OA PRG4+ lining fibroblasts"

oxphos_markers <- c("NDUFS3", "NDUFS1", "NDUFS2", "NDUFS7", "NDUFS8", "NDUFV1", "NDUFV2", "NDUFA9", "ATP5F", "COX5A", "UQCRC1")

min_cells_per_subtype <- 100

min_cells_per_gene <- 10

max_genes_use <- 1500

max_cells_use <- 3000

nc_nNet_use <- 10

nc_nComp_use <- 3

nc_q_use <- 0.9

td_K_use <- 3

ma_nDim_use <- 2

n_cores_use <- 2L

Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "2")

set.seed(1234)

message2 <- function(...) {
    message(paste0(...))
}

save_plot_pdf_png <- function(plot, pdf_file, png_file, width, height, dpi = 300) {
    ggsave(filename = pdf_file, plot = plot, width = width, height = height, limitsize = FALSE)
    ggsave(filename = png_file, plot = plot, width = width, height = height, dpi = dpi, limitsize = FALSE)
}

get_counts_matrix <- function(obj, assay = "RNA", target_cells = NULL) {
    if (!assay %in% SeuratObject::Assays(obj)) {
        stop("对象中不存在 assay：", assay)
    }
    assay_obj <- obj[[assay]]
    layer_names <- SeuratObject::Layers(assay_obj)
    count_layers <- layer_names[grepl("^counts($|\\.)", layer_names)]
    if (length(count_layers) == 0L) {
        stop("assay ", assay, " 中没有 counts layer。")
    }
    mats <- lapply(count_layers, function(layer_name) {
        x <- SeuratObject::LayerData(assay_obj, layer = layer_name)
        if (!is.null(target_cells)) {
            keep <- intersect(target_cells, colnames(x))
            if (length(keep) == 0L) 
                return(NULL)
            x <- x[, keep, drop = FALSE]
        }
        x
    })
    mats <- Filter(Negate(is.null), mats)
    if (length(mats) == 0L) 
        stop("counts layer 中没有匹配的细胞。")
    common_genes <- Reduce(intersect, lapply(mats, rownames))
    mats <- lapply(mats, function(x) x[common_genes, , drop = FALSE])
    counts <- do.call(cbind, mats)
    if (anyDuplicated(colnames(counts))) {
        stop("counts 矩阵中存在重复条形码。")
    }
    if (!is.null(target_cells)) {
        missing_cells <- setdiff(target_cells, colnames(counts))
        if (length(missing_cells) > 0L) {
            stop("counts 矩阵缺少 ", length(missing_cells), " 个目标细胞。")
        }
        counts <- counts[, target_cells, drop = FALSE]
    }
    methods::as(counts, "dgCMatrix")
}

clean_name <- function(x) {
    gsub("[^A-Za-z0-9_]+", "_", x)
}

prepare_knk_matrix <- function(full_counts_mat, full_meta, cells_use, subtype_name, target_gene = "NDUFS3", oxphos_markers = NULL, 
    max_genes_use = 3000, max_cells_use = 3000, min_cells_per_gene = 10) {
    message2("开始准备 ", subtype_name, " 的 scTenifoldKnk 输入矩阵。")
    cells_use <- intersect(cells_use, colnames(full_counts_mat))
    if (length(cells_use) == 0) {
        stop(subtype_name, " 在 counts 矩阵中没有匹配到细胞。")
    }
    counts_sub <- full_counts_mat[, cells_use, drop = FALSE]
    meta_sub <- full_meta[cells_use, , drop = FALSE]
    stopifnot(identical(colnames(counts_sub), rownames(meta_sub)))
    if (is.finite(max_cells_use) && ncol(counts_sub) > max_cells_use) {
        set.seed(1234)
        sampled_cells <- sample(colnames(counts_sub), size = max_cells_use, replace = FALSE)
        counts_sub <- counts_sub[, sampled_cells, drop = FALSE]
        meta_sub <- meta_sub[sampled_cells, , drop = FALSE]
        message2(subtype_name, " 细胞数超过 max_cells_use，已随机抽样到：", ncol(counts_sub))
    }
    if (!target_gene %in% rownames(counts_sub)) {
        stop(subtype_name, " 中没有找到目标基因：", target_gene)
    }
    seu_obj <- CreateSeuratObject(counts = counts_sub, meta.data = meta_sub, assay = "RNA")
    DefaultAssay(seu_obj) <- "RNA"
    seu_obj <- NormalizeData(seu_obj, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
    seu_obj <- FindVariableFeatures(seu_obj, selection.method = "vst", nfeatures = max_genes_use, verbose = FALSE)
    variable_genes <- VariableFeatures(seu_obj)
    counts_mat <- get_counts_matrix(seu_obj, assay = "RNA")
    gene_detected_cells <- Matrix::rowSums(counts_mat > 0)
    expressed_genes <- names(gene_detected_cells)[gene_detected_cells >= min_cells_per_gene]
    genes_use <- unique(c(variable_genes, target_gene, oxphos_markers))
    genes_use <- genes_use[genes_use %in% rownames(counts_mat)]
    genes_use <- genes_use[genes_use %in% expressed_genes | genes_use == target_gene]
    if (!target_gene %in% genes_use) {
        genes_use <- unique(c(target_gene, genes_use))
    }
    if (length(genes_use) > max_genes_use) {
        keep_priority <- unique(c(target_gene, oxphos_markers[oxphos_markers %in% genes_use]))
        keep_other <- setdiff(genes_use, keep_priority)
        n_other_keep <- max(0, max_genes_use - length(keep_priority))
        genes_use <- unique(c(keep_priority, head(keep_other, n_other_keep)))
    }
    genes_use <- intersect(genes_use, rownames(counts_mat))
    counts_mat_use <- counts_mat[genes_use, , drop = FALSE]
    ndufs3_positive_cells <- sum(counts_mat_use[target_gene, ] > 0)
    message2(subtype_name, " 中 ", target_gene, " 阳性细胞数：", ndufs3_positive_cells, " / ", ncol(counts_mat_use))
    if (ndufs3_positive_cells == 0) {
        stop(subtype_name, " 中 NDUFS3 表达为 0，无法进行可靠的 NDUFS3 虚拟敲除。")
    }
    counts_mat_use <- as.matrix(counts_mat_use)
    storage.mode(counts_mat_use) <- "numeric"
    return(list(seu_obj = seu_obj, count_matrix = counts_mat_use, genes_use = genes_use, ndufs3_positive_cells = ndufs3_positive_cells))
}

extract_diff_regulation <- function(knk_res, subtype_name, target_gene = "NDUFS3") {
    if (!"diffRegulation" %in% names(knk_res)) {
        stop("scTenifoldKnk 结果中没有 diffRegulation。")
    }
    diff_df <- knk_res$diffRegulation %>% as.data.frame()
    if (!"gene" %in% colnames(diff_df)) {
        diff_df <- diff_df %>% tibble::rownames_to_column("gene")
    }
    if (!"p.adj" %in% colnames(diff_df)) {
        if ("p_adj" %in% colnames(diff_df)) {
            diff_df$p.adj <- diff_df$p_adj
        }
        else if ("padj" %in% colnames(diff_df)) {
            diff_df$p.adj <- diff_df$padj
        }
        else {
            diff_df$p.adj <- NA
        }
    }
    if (!"p.value" %in% colnames(diff_df)) {
        if ("p_val" %in% colnames(diff_df)) {
            diff_df$p.value <- diff_df$p_val
        }
        else {
            diff_df$p.value <- NA
        }
    }
    if (!"Z" %in% colnames(diff_df)) {
        diff_df$Z <- NA
    }
    if (!"distance" %in% colnames(diff_df)) {
        diff_df$distance <- NA
    }
    if (!"FC" %in% colnames(diff_df)) {
        diff_df$FC <- NA
    }
    diff_df <- diff_df %>% dplyr::mutate(subtype = subtype_name, target_gene = target_gene, neg_log10_padj = -log10(p.adj + 
        9.99999998481683e-301), significant = ifelse(!is.na(p.adj) & p.adj < 0.05, "Significant", "Not significant")) %>% 
        dplyr::arrange(p.adj, dplyr::desc(abs(Z)))
    diff_df$rank <- seq_len(nrow(diff_df))
    return(diff_df)
}

plot_ko_network <- function(knk_res, subtype_name, target_gene = "NDUFS3", out_prefix) {
    message2("开始绘制 ", subtype_name, " 的 KO network 图。")
    pdf_file <- paste0("05_KO_network_plot/", out_prefix, "_NDUFS3_KO_network.pdf")
    png_file <- paste0("05_KO_network_plot/", out_prefix, "_NDUFS3_KO_network.png")
    tryCatch({
        pdf(pdf_file, width = 8, height = 8)
        plotKO(X = knk_res, gKO = target_gene, q = 0.99, annotate = FALSE, fdrThreshold = 0.05)
        dev.off()
    }, error = function(e) {
        if (dev.cur() != 1) 
            dev.off()
        warning("PDF KO network 绘制失败：", e$message)
    })
    tryCatch({
        png(png_file, width = 2400, height = 2400, res = 300)
        plotKO(X = knk_res, gKO = target_gene, q = 0.99, annotate = FALSE, fdrThreshold = 0.05)
        dev.off()
    }, error = function(e) {
        if (dev.cur() != 1) 
            dev.off()
        warning("PNG KO network 绘制失败：", e$message)
    })
}

plot_volcano <- function(diff_df, subtype_name, out_prefix, oxphos_markers = NULL, top_n_label = 15) {
    label_genes <- diff_df %>% dplyr::filter(!is.na(p.adj)) %>% dplyr::arrange(p.adj, dplyr::desc(abs(Z))) %>% dplyr::slice_head(n = top_n_label) %>% 
        dplyr::pull(gene)
    label_genes <- unique(c(label_genes, oxphos_markers[oxphos_markers %in% diff_df$gene], "NDUFS3"))
    p <- ggplot(diff_df, aes(x = Z, y = neg_log10_padj)) + geom_point(aes(color = significant), alpha = 0.75, size = 1.4) + 
        geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.3) + geom_hline(yintercept = -log10(0.05), linetype = "dashed", 
        linewidth = 0.3) + ggrepel::geom_text_repel(data = diff_df %>% dplyr::filter(gene %in% label_genes), aes(label = gene), 
        size = 3, max.overlaps = 50, box.padding = 0.35, point.padding = 0.25) + scale_color_manual(values = c(Significant = "#D73027", 
        `Not significant` = "grey70")) + theme_bw(base_size = 13) + labs(title = paste0(subtype_name, ": NDUFS3 virtual knockout"), 
        x = "Differential regulation Z-score", y = "-log10(adjusted p value)", color = NULL) + theme(plot.title = element_text(hjust = 0.5, 
        face = "bold"), panel.grid = element_blank())
    save_plot_pdf_png(p, paste0("06_volcano_plot/", out_prefix, "_NDUFS3_KO_volcano.pdf"), paste0("06_volcano_plot/", out_prefix, 
        "_NDUFS3_KO_volcano.png"), width = 7, height = 6)
    return(p)
}

plot_top_genes <- function(diff_df, subtype_name, out_prefix, top_n = 30) {
    top_df <- diff_df %>% dplyr::filter(!is.na(p.adj)) %>% dplyr::arrange(p.adj, dplyr::desc(abs(Z))) %>% dplyr::slice_head(n = top_n) %>% 
        dplyr::mutate(gene = factor(gene, levels = rev(gene)))
    p <- ggplot(top_df, aes(x = gene, y = Z)) + geom_col(width = 0.75) + coord_flip() + theme_bw(base_size = 13) + labs(title = paste0(subtype_name, 
        ": Top ", top_n, " perturbed genes after NDUFS3 KO"), x = NULL, y = "Differential regulation Z-score") + theme(plot.title = element_text(hjust = 0.5, 
        face = "bold"), panel.grid = element_blank(), axis.text.y = element_text(face = "italic"))
    save_plot_pdf_png(p, paste0("07_top_gene_plot/", out_prefix, "_NDUFS3_KO_top", top_n, "_genes.pdf"), paste0("07_top_gene_plot/", 
        out_prefix, "_NDUFS3_KO_top", top_n, "_genes.png"), width = 7, height = 8)
    return(p)
}

prepare_manifold_df <- function(knk_res, diff_df) {
    if (!"manifoldAlignment" %in% names(knk_res)) {
        warning("scTenifoldKnk 结果中没有 manifoldAlignment。")
        return(NULL)
    }
    ma <- knk_res$manifoldAlignment %>% as.data.frame()
    if (ncol(ma) < 2) {
        warning("manifoldAlignment 维度少于 2，无法绘图。")
        return(NULL)
    }
    colnames(ma)[1:2] <- c("Dim1", "Dim2")
    ma$row_id <- rownames(ma)
    n_gene <- nrow(diff_df)
    if (all(grepl("^[XY]_", ma$row_id))) {
        ma$gene <- sub("^[XY]_", "", ma$row_id)
        ma$condition <- ifelse(grepl("^Y_", ma$row_id), "NDUFS3_KO", "Original")
    }
    else {
        ma$condition <- dplyr::case_when(grepl("KO|knock|Y", ma$row_id, ignore.case = TRUE) ~ "NDUFS3_KO", TRUE ~ "Original")
        ma$gene <- ma$row_id %>% gsub("_KO|\\.KO|-KO", "", ., ignore.case = TRUE) %>% gsub("_WT|\\.WT|-WT", "", ., ignore.case = TRUE) %>% 
            gsub("_X|\\.X|-X", "", ., ignore.case = TRUE) %>% gsub("_Y|\\.Y|-Y", "", ., ignore.case = TRUE)
    }
    ma <- ma %>% dplyr::left_join(diff_df %>% dplyr::select(gene, Z, distance, p.value, p.adj, significant), by = "gene")
    return(ma)
}

plot_manifold <- function(knk_res, diff_df, subtype_name, out_prefix, top_n_label = 20) {
    ma <- prepare_manifold_df(knk_res, diff_df)
    if (is.null(ma)) {
        return(NULL)
    }
    top_genes <- diff_df %>% dplyr::filter(!is.na(p.adj)) %>% dplyr::arrange(p.adj, dplyr::desc(abs(Z))) %>% dplyr::slice_head(n = top_n_label) %>% 
        dplyr::pull(gene)
    ma_top <- ma %>% dplyr::filter(gene %in% top_genes)
    p <- ggplot(ma, aes(x = Dim1, y = Dim2)) + geom_point(aes(color = condition), alpha = 0.45, size = 1.2) + ggrepel::geom_text_repel(data = ma_top, 
        aes(label = gene), size = 3, max.overlaps = 60, box.padding = 0.35) + theme_bw(base_size = 13) + labs(title = paste0(subtype_name, 
        ": manifold alignment after NDUFS3 KO"), x = "Manifold dimension 1", y = "Manifold dimension 2", color = NULL) + 
        theme(plot.title = element_text(hjust = 0.5, face = "bold"), panel.grid = element_blank())
    save_plot_pdf_png(p, paste0("08_manifold_plot/", out_prefix, "_NDUFS3_KO_manifold_alignment.pdf"), paste0("08_manifold_plot/", 
        out_prefix, "_NDUFS3_KO_manifold_alignment.png"), width = 7, height = 6)
    write.csv(ma, file = paste0("08_manifold_plot/", out_prefix, "_NDUFS3_KO_manifold_alignment_coordinates.csv"), row.names = FALSE)
    return(p)
}

plot_marker_response <- function(diff_df, subtype_name, out_prefix, markers) {
    marker_df <- diff_df %>% dplyr::filter(gene %in% markers)
    if (nrow(marker_df) == 0) {
        warning(subtype_name, " 中没有找到指定炎症 marker 的 diffRegulation 结果。")
        return(NULL)
    }
    marker_levels <- rev(markers[markers %in% marker_df$gene])
    marker_df <- marker_df %>% dplyr::mutate(gene = factor(gene, levels = marker_levels), label = ifelse(is.na(p.adj), "NA", 
        paste0("FDR=", signif(p.adj, 3))))
    write.csv(marker_df, file = paste0("09_marker_response_plot/", out_prefix, "_NDUFS3_KO_oxphos_marker_response_table.csv"), 
        row.names = FALSE)
    p <- ggplot(marker_df, aes(x = gene, y = Z)) + geom_col(width = 0.7) + geom_text(aes(label = label), hjust = ifelse(marker_df$Z >= 
        0, -0.05, 1.05), size = 3) + coord_flip() + theme_bw(base_size = 13) + labs(title = paste0(subtype_name, ": OXPHOS/complex I marker response after NDUFS3 KO"), 
        x = NULL, y = "Differential regulation Z-score") + theme(plot.title = element_text(hjust = 0.5, face = "bold"), panel.grid = element_blank(), 
        axis.text.y = element_text(face = "italic"))
    save_plot_pdf_png(p, paste0("09_marker_response_plot/", out_prefix, "_NDUFS3_KO_oxphos_marker_response.pdf"), paste0("09_marker_response_plot/", 
        out_prefix, "_NDUFS3_KO_oxphos_marker_response.png"), width = 7, height = 5.8)
    return(p)
}

run_go_enrichment <- function(diff_df, subtype_name, out_prefix, top_n_for_enrich = 200) {
    sig_genes <- diff_df %>% dplyr::filter(!is.na(p.adj), p.adj < 0.05) %>% dplyr::pull(gene) %>% unique()
    if (length(sig_genes) < 10) {
        warning(subtype_name, " 显著扰动基因少于 10 个，改用 Top ", top_n_for_enrich, " 基因做探索性富集。")
        sig_genes <- diff_df %>% dplyr::filter(!is.na(p.adj)) %>% dplyr::arrange(p.adj, dplyr::desc(abs(Z))) %>% dplyr::slice_head(n = top_n_for_enrich) %>% 
            dplyr::pull(gene) %>% unique()
    }
    sig_genes <- sig_genes[!is.na(sig_genes)]
    if (length(sig_genes) < 5) {
        warning(subtype_name, " 可用于富集分析的基因少于 5 个，跳过 GO 富集。")
        return(NULL)
    }
    gene_map <- AnnotationDbi::select(org.Hs.eg.db, keys = sig_genes, columns = c("SYMBOL", "ENTREZID"), keytype = "SYMBOL") %>% 
        dplyr::filter(!is.na(ENTREZID)) %>% dplyr::distinct(SYMBOL, ENTREZID)
    if (nrow(gene_map) < 5) {
        warning(subtype_name, " SYMBOL 到 ENTREZID 映射后少于 5 个基因，跳过 GO 富集。")
        return(NULL)
    }
    ego <- enrichGO(gene = gene_map$ENTREZID, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "BP", pAdjustMethod = "BH", 
        pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = TRUE)
    ego_df <- as.data.frame(ego)
    write.csv(ego_df, file = paste0("10_GO_enrichment/", out_prefix, "_NDUFS3_KO_GO_BP_enrichment.csv"), row.names = FALSE)
    if (nrow(ego_df) == 0) {
        warning(subtype_name, " 没有显著 GO BP 富集结果。")
        return(NULL)
    }
    p <- enrichplot::dotplot(ego, showCategory = min(20, nrow(ego_df))) + ggtitle(paste0(subtype_name, ": GO BP enrichment of NDUFS3-KO perturbed genes")) + 
        theme_bw(base_size = 12) + theme(plot.title = element_text(hjust = 0.5, face = "bold"))
    save_plot_pdf_png(p, paste0("10_GO_enrichment/", out_prefix, "_NDUFS3_KO_GO_BP_dotplot.pdf"), paste0("10_GO_enrichment/", 
        out_prefix, "_NDUFS3_KO_GO_BP_dotplot.png"), width = 9, height = 7)
    return(list(ego = ego, ego_df = ego_df, gene_map = gene_map))
}

run_ndufs3_knockout_for_subtype <- function(sc_fibroblast, full_counts_mat, full_meta, subtype_name, target_gene = "NDUFS3", 
    oxphos_markers = oxphos_markers, min_cells_per_subtype = 100, min_cells_per_gene = 10, max_genes_use = 3000, max_cells_use = 3000, 
    nc_nNet_use = 10, nc_nComp_use = 3, nc_q_use = 0.9, td_K_use = 3, ma_nDim_use = 2, n_cores_use = 4) {
    message2("\n====================================================")
    message2("开始分析亚群：", subtype_name)
    message2("====================================================\n")
    out_prefix <- clean_name(subtype_name)
    cells_use <- rownames(full_meta)[as.character(full_meta$Group) == "OA" & as.character(full_meta$fibroblast_subtype_published) == 
        "PRG4+ lining fibroblasts"]
    cells_use <- intersect(cells_use, colnames(full_counts_mat))
    if (length(cells_use) < min_cells_per_subtype) {
        warning(subtype_name, " 细胞数少于 ", min_cells_per_subtype, "，跳过该亚群。")
        return(NULL)
    }
    prep <- prepare_knk_matrix(full_counts_mat = full_counts_mat, full_meta = full_meta, cells_use = cells_use, subtype_name = subtype_name, 
        target_gene = target_gene, oxphos_markers = oxphos_markers, max_genes_use = max_genes_use, max_cells_use = max_cells_use, 
        min_cells_per_gene = min_cells_per_gene)
    sc_sub_clean <- prep$seu_obj
    saveRDS(sc_sub_clean, file = paste0("01_object/sc_fibroblast_", out_prefix, "_clean_for_NDUFS3_KO.rds"))
    save(sc_sub_clean, file = paste0("01_object/sc_fibroblast_", out_prefix, "_clean_for_NDUFS3_KO.rda"))
    count_matrix <- prep$count_matrix
    message2(subtype_name, " 输入矩阵维度：")
    print(dim(count_matrix))
    saveRDS(count_matrix, file = paste0("02_count_matrix/", out_prefix, "_NDUFS3_KO_input_count_matrix.rds"))
    write.table(rownames(count_matrix), file = paste0("02_count_matrix/", out_prefix, "_NDUFS3_KO_input_genes.txt"), quote = FALSE, 
        row.names = FALSE, col.names = FALSE)
    n_cells_now <- ncol(count_matrix)
    nc_nCells_use <- min(500, max(50, floor(n_cells_now * 0.8)))
    if (n_cells_now < 80) {
        nc_nCells_use <- max(20, floor(n_cells_now * 0.8))
    }
    message2(subtype_name, " scTenifoldKnk nc_nCells_use = ", nc_nCells_use)
    model_rds <- paste0("03_scTenifoldKnk_result/", out_prefix, "_NDUFS3_KO_scTenifoldKnk_result.rds")
    if (file.exists(model_rds)) {
        message2("检测到已完成模型，直接读取：", model_rds)
        knk_res <- readRDS(model_rds)
    }
    else {
        knk_res <- scTenifoldKnk(countMatrix = count_matrix, qc = FALSE, gKO = target_gene, qc_mtThreshold = 0.2, qc_minLSize = 0, 
            qc_minCells = min_cells_per_gene, nc_lambda = 0, nc_nNet = nc_nNet_use, nc_nCells = nc_nCells_use, nc_nComp = nc_nComp_use, 
            nc_scaleScores = TRUE, nc_symmetric = FALSE, nc_q = nc_q_use, td_K = td_K_use, td_maxIter = 1000, td_maxError = 1e-05, 
            td_nDecimal = 3, ma_nDim = ma_nDim_use, nCores = n_cores_use)
        saveRDS(knk_res, model_rds, compress = FALSE)
        save(knk_res, file = paste0("03_scTenifoldKnk_result/", out_prefix, "_NDUFS3_KO_scTenifoldKnk_result.rda"))
    }
    diff_df <- extract_diff_regulation(knk_res = knk_res, subtype_name = subtype_name, target_gene = target_gene)
    write.csv(diff_df, file = paste0("04_diffRegulation_table/", out_prefix, "_NDUFS3_KO_diffRegulation_all_genes.csv"), 
        row.names = FALSE)
    sig_df <- diff_df %>% dplyr::filter(!is.na(p.adj), p.adj < 0.05)
    write.csv(sig_df, file = paste0("04_diffRegulation_table/", out_prefix, "_NDUFS3_KO_diffRegulation_sig_FDR0.05.csv"), 
        row.names = FALSE)
    top30_df <- diff_df %>% dplyr::filter(!is.na(p.adj)) %>% dplyr::arrange(p.adj, dplyr::desc(abs(Z))) %>% dplyr::slice_head(n = 30)
    write.csv(top30_df, file = paste0("04_diffRegulation_table/", out_prefix, "_NDUFS3_KO_top30_perturbed_genes.csv"), row.names = FALSE)
    plot_ko_network(knk_res = knk_res, subtype_name = subtype_name, target_gene = target_gene, out_prefix = out_prefix)
    p_volcano <- plot_volcano(diff_df = diff_df, subtype_name = subtype_name, out_prefix = out_prefix, oxphos_markers = oxphos_markers, 
        top_n_label = 15)
    p_top <- plot_top_genes(diff_df = diff_df, subtype_name = subtype_name, out_prefix = out_prefix, top_n = 30)
    p_manifold <- plot_manifold(knk_res = knk_res, diff_df = diff_df, subtype_name = subtype_name, out_prefix = out_prefix, 
        top_n_label = 20)
    p_marker <- plot_marker_response(diff_df = diff_df, subtype_name = subtype_name, out_prefix = out_prefix, markers = oxphos_markers)
    go_res <- run_go_enrichment(diff_df = diff_df, subtype_name = subtype_name, out_prefix = out_prefix, top_n_for_enrich = 200)
    summary_df <- data.frame(subtype = subtype_name, target_gene = target_gene, n_cells_original = length(cells_use), n_cells_used = ncol(count_matrix), 
        n_genes_used = nrow(count_matrix), ndufs3_positive_cells = prep$ndufs3_positive_cells, ndufs3_positive_ratio = prep$ndufs3_positive_cells/ncol(count_matrix), 
        n_sig_genes_FDR0.05 = nrow(sig_df), top_gene = ifelse(nrow(diff_df) > 0, diff_df$gene[1], NA), top_gene_Z = ifelse(nrow(diff_df) > 
            0, diff_df$Z[1], NA), top_gene_padj = ifelse(nrow(diff_df) > 0, diff_df$p.adj[1], NA))
    write.csv(summary_df, file = paste0("11_summary_tables/", out_prefix, "_NDUFS3_KO_summary.csv"), row.names = FALSE)
    message2("完成亚群：", subtype_name)
    return(list(subtype = subtype_name, sc_sub_clean = sc_sub_clean, count_matrix = count_matrix, knk_res = knk_res, diff_df = diff_df, 
        sig_df = sig_df, top30_df = top30_df, summary_df = summary_df, go_res = go_res))
}

if (!file.exists(input_rds)) {
    stop("没有找到最终成纤维细胞 RDS：", input_rds)
}

if (!file.exists(whole_atlas_rds)) {
    stop("没有找到全图谱 RDS：", whole_atlas_rds)
}

sc_fibroblast <- readRDS(input_rds)

whole_atlas <- readRDS(whole_atlas_rds)

if (!inherits(sc_fibroblast, "Seurat")) {
    stop("最终成纤维细胞输入不是 Seurat 对象。")
}

if (!inherits(whole_atlas, "Seurat")) {
    stop("全图谱输入不是 Seurat 对象。")
}

if (!"fibroblast_subtype_published" %in% colnames(sc_fibroblast@meta.data)) {
    stop("最终对象缺少 fibroblast_subtype_published 注释列。")
}

need_cols <- c("fibroblast_subtype_published", "Group", "orig.ident")

missing_cols <- setdiff(need_cols, colnames(sc_fibroblast@meta.data))

if (length(missing_cols) > 0) {
    stop("最终成纤维细胞对象缺少以下 metadata 列：", paste(missing_cols, collapse = ", "))
}

sc_fibroblast$fibroblast_subtype_published <- as.character(sc_fibroblast$fibroblast_subtype_published)

sc_fibroblast$fibroblast_subtype_published <- factor(sc_fibroblast$fibroblast_subtype_published, levels = locked_fibroblast_subtypes)

sc_fibroblast$Group <- factor(as.character(sc_fibroblast$Group), levels = c("Control", "OA"))

if (anyNA(sc_fibroblast$fibroblast_subtype_published)) {
    stop("存在不属于最终7个锁定亚群的细胞。")
}

if (anyNA(sc_fibroblast$Group)) {
    stop("Group 中存在非 Control/OA 标签。")
}

cat("最终成纤维细胞亚群数量：\n")

print(table(sc_fibroblast$fibroblast_subtype_published, useNA = "ifany"))

cat("最终成纤维细胞亚群与 Group：\n")

print(table(sc_fibroblast$fibroblast_subtype_published, sc_fibroblast$Group, useNA = "ifany"))

cat("样本来源：\n")

print(table(sc_fibroblast$orig.ident, useNA = "ifany"))

subtype_count_df <- sc_fibroblast@meta.data %>% dplyr::count(fibroblast_subtype_published, Group)

write.csv(subtype_count_df, file = "11_summary_tables/input_Fibroblast_subtype_by_Group_cell_counts.csv", row.names = FALSE)

fibro_cells <- colnames(sc_fibroblast)

full_meta <- sc_fibroblast@meta.data[fibro_cells, , drop = FALSE]

full_counts_mat <- NULL

count_source <- NULL

if ("RNA" %in% SeuratObject::Assays(whole_atlas)) {
    full_counts_mat <- tryCatch(get_counts_matrix(whole_atlas, assay = "RNA", target_cells = fibro_cells), error = function(e) {
        message("全图谱 RNA counts 不可用：", e$message)
        NULL
    })
    if (!is.null(full_counts_mat)) {
        count_source <- "whole-atlas RNA counts"
    }
}

if (is.null(full_counts_mat)) {
    full_counts_mat <- get_counts_matrix(sc_fibroblast, assay = "SCT", target_cells = fibro_cells)
    count_source <- "final-fibroblast SCT counts layer"
}

if (!target_gene %in% rownames(full_counts_mat)) {
    stop("实际用于虚拟敲除的 counts 矩阵中没有找到目标基因：", target_gene)
}

if (!identical(colnames(full_counts_mat), rownames(full_meta))) {
    full_counts_mat <- full_counts_mat[, rownames(full_meta), drop = FALSE]
}

if (!identical(colnames(full_counts_mat), rownames(full_meta))) {
    stop("counts 矩阵与 metadata 条形码顺序不一致。")
}

write.csv(data.frame(metric = c("fibroblast_cells", "count_genes", "count_source"), value = c(ncol(full_counts_mat), nrow(full_counts_mat), 
    count_source)), file = "11_summary_tables/input_count_source.csv", row.names = FALSE)

message("完整成纤维细胞 counts 矩阵维度：")

print(dim(full_counts_mat))

message("counts 来源：", count_source)

rm(whole_atlas)

invisible(gc())

all_results <- list()

for (subtype_now in target_subtypes) {
    res_now <- run_ndufs3_knockout_for_subtype(sc_fibroblast = sc_fibroblast, full_counts_mat = full_counts_mat, full_meta = full_meta, 
        subtype_name = subtype_now, target_gene = target_gene, oxphos_markers = oxphos_markers, min_cells_per_subtype = min_cells_per_subtype, 
        min_cells_per_gene = min_cells_per_gene, max_genes_use = max_genes_use, max_cells_use = max_cells_use, nc_nNet_use = nc_nNet_use, 
        nc_nComp_use = nc_nComp_use, nc_q_use = nc_q_use, td_K_use = td_K_use, ma_nDim_use = ma_nDim_use, n_cores_use = n_cores_use)
    all_results[[subtype_now]] <- res_now
}

summary_all <- purrr::map_dfr(all_results, function(x) {
    if (is.null(x)) {
        return(NULL)
    }
    x$summary_df
})

write.csv(summary_all, file = "11_summary_tables/NDUFS3_KO_all_subtypes_summary.csv", row.names = FALSE)

diff_all <- purrr::map_dfr(all_results, function(x) {
    if (is.null(x)) {
        return(NULL)
    }
    x$diff_df
})

write.csv(diff_all, file = "11_summary_tables/NDUFS3_KO_all_subtypes_diffRegulation_merged.csv", row.names = FALSE)

sig_all <- diff_all %>% dplyr::filter(!is.na(p.adj), p.adj < 0.05)

write.csv(sig_all, file = "11_summary_tables/NDUFS3_KO_all_subtypes_sig_genes_FDR0.05_merged.csv", row.names = FALSE)

top_compare <- diff_all %>% dplyr::filter(!is.na(p.adj)) %>% dplyr::group_by(subtype) %>% dplyr::arrange(p.adj, dplyr::desc(abs(Z)), 
    .by_group = TRUE) %>% dplyr::slice_head(n = 20) %>% dplyr::ungroup()

write.csv(top_compare, file = "11_summary_tables/NDUFS3_KO_top20_genes_each_subtype.csv", row.names = FALSE)

if (nrow(top_compare) > 0) {
    p_top_compare <- ggplot(top_compare, aes(x = reorder(gene, Z), y = Z)) + geom_col(width = 0.75) + coord_flip() + facet_wrap(~subtype, 
        scales = "free_y") + theme_bw(base_size = 13) + labs(title = "Top perturbed genes after NDUFS3 virtual knockout", 
        x = NULL, y = "Differential regulation Z-score") + theme(plot.title = element_text(hjust = 0.5, face = "bold"), panel.grid = element_blank(), 
        axis.text.y = element_text(face = "italic"), strip.text = element_text(face = "bold"))
    save_plot_pdf_png(p_top_compare, "11_summary_tables/NDUFS3_KO_top20_genes_all_subtypes_compare.pdf", "11_summary_tables/NDUFS3_KO_top20_genes_all_subtypes_compare.png", 
        width = 12, height = 8)
}

marker_compare <- diff_all %>% dplyr::filter(gene %in% oxphos_markers) %>% dplyr::mutate(gene = factor(gene, levels = oxphos_markers))

write.csv(marker_compare, file = "11_summary_tables/NDUFS3_KO_oxphos_marker_response_all_subtypes.csv", row.names = FALSE)

if (nrow(marker_compare) > 0) {
    p_marker_compare <- ggplot(marker_compare, aes(x = gene, y = Z, fill = subtype)) + geom_col(position = position_dodge(width = 0.75), 
        width = 0.65) + geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3) + theme_bw(base_size = 13) + labs(title = "OXPHOS/complex I marker response after NDUFS3 virtual knockout", 
        x = NULL, y = "Differential regulation Z-score", fill = NULL) + theme(plot.title = element_text(hjust = 0.5, face = "bold"), 
        panel.grid = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"))
    save_plot_pdf_png(p_marker_compare, "11_summary_tables/NDUFS3_KO_oxphos_marker_response_all_subtypes_compare.pdf", "11_summary_tables/NDUFS3_KO_oxphos_marker_response_all_subtypes_compare.png", 
        width = 9, height = 5.5)
}

saveRDS(all_results, file = "11_summary_tables/NDUFS3_KO_all_results_list.rds")

save(all_results, summary_all, diff_all, sig_all, top_compare, marker_compare, file = "11_summary_tables/NDUFS3_KO_all_results_workspace.rda")

message("NDUFS3 在OA组 PRG4+ lining fibroblasts 中的虚拟敲除分析全部完成。")

