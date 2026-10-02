args <- commandArgs(trailingOnly = FALSE)

script_file <- sub("^--file=", "", grep("^--file=", args, value = TRUE)[1])

code_root <- dirname(dirname(normalizePath(script_file, mustWork = TRUE)))

source(file.path(code_root, "00_configuration.R"))

local({
    INPUT_DIR <- file.path(input_root, "singlecell")
    OUTPUT_DIR <- file.path(output_root, "39_NDUFS3_fibroblast_Monocle2")
    FINAL_FIBRO_RDS <- file.path(INPUT_DIR, "GSE216651_Fibroblasts_final_7_published_subtypes.rds")
    WHOLE_ATLAS_RDS <- file.path(INPUT_DIR, "GSE216651_scRNA_PC30_res0.6_without_clusters16_17.rds")
    OXPHOS_SCORE_CSV <- paste0(paste0(file.path(output_root, "33_phosphorylation_score_major_and_fibroblast"), "/"), "02_final_fibroblast_phosphorylation_score_5methods.csv")
    RDS_DIR <- file.path(OUTPUT_DIR, "01_rds")
    TABLE_DIR <- file.path(OUTPUT_DIR, "02_tables")
    FIGURE_DIR <- file.path(OUTPUT_DIR, "03_figures_PDF_PNG")
    LOG_DIR <- file.path(OUTPUT_DIR, "04_logs")
    for (x in c(RDS_DIR, TABLE_DIR, FIGURE_DIR, LOG_DIR)) {
        dir.create(x, recursive = TRUE, showWarnings = FALSE)
    }
    GROUP_FIELD <- "Group"
    SUBTYPE_FIELD <- "fibroblast_subtype_published"
    TARGET_GENE <- "NDUFS3"
    PRG4_SUBTYPE <- "PRG4+ lining fibroblasts"
    EXPECTED_MONOCLE_VERSION <- "2.36.9213"
    EXPECTED_FIBRO_CELLS <- 31004L
    EXPECTED_ATLAS_CELLS <- 60302L
    CONDITION_LEVELS <- c("Control", "OA")
    SUBTYPE_LEVELS <- c("PRG4+ lining fibroblasts", "APOD+ fibroblasts", "CD34+ sublining fibroblasts", "CXCL12+ sublining fibroblasts", 
        "DKK3+ sublining fibroblasts", "RSPO3+ fibroblasts", "POSTN+ fibroblasts")
    EXPECTED_SUBTYPE_COUNTS <- c(`PRG4+ lining fibroblasts` = 9190L, `APOD+ fibroblasts` = 6161L, `CD34+ sublining fibroblasts` = 4714L, 
        `CXCL12+ sublining fibroblasts` = 5357L, `DKK3+ sublining fibroblasts` = 3220L, `RSPO3+ fibroblasts` = 1749L, `POSTN+ fibroblasts` = 613L)
    SUBTYPE_COLORS <- c(`PRG4+ lining fibroblasts` = "#4E79A7", `APOD+ fibroblasts` = "#59A14F", `CD34+ sublining fibroblasts` = "#F28E2B", 
        `CXCL12+ sublining fibroblasts` = "#E15759", `DKK3+ sublining fibroblasts` = "#B07AA1", `RSPO3+ fibroblasts` = "#FF9DA7", 
        `POSTN+ fibroblasts` = "#9C755F")
    GROUP_COLORS <- c(Control = "#2878B5", OA = "#E3B341")
    NDUFS3_COLORS <- c(`NDUFS3+` = "#E64B35", `NDUFS3-` = "#4DBBD5")
    PSEUDOTIME_COLORS <- c("#482878", "#31688E", "#35B779", "#FDE725")
    OXPHOS_COLORS <- c("#440154", "#31688E", "#35B779", "#FDE725")
    MIN_PCT <- 0.1
    LOGFC_THRESHOLD <- 0.25
    ORDERING_FDR <- 0.05
    MAX_ORDERING_GENES <- 1500L
    MIN_ORDERING_GENES <- 50L
    TOP_HEATMAP_GENES <- 50L
    CORES <- 1L
    SEED <- 20260807L
    REUSE_MARKER_TABLES <- TRUE
    REUSE_SAVED_CDS <- TRUE
    log_file <- file.path(LOG_DIR, "39_NDUFS3_fibroblast_Monocle2_run.log")
    log_con <- file(log_file, open = "wt", encoding = "UTF-8")
    sink(log_con, type = "output", split = TRUE)
    on.exit({
        try(sink(type = "output"), silent = TRUE)
        try(close(log_con), silent = TRUE)
    }, add = TRUE)
    options(stringsAsFactors = FALSE)
    options(future.globals.maxSize = 32 * 1024^3)
    set.seed(SEED)
    message("Analysis started: ", format(Sys.time()))
    message("Output: ", OUTPUT_DIR)
    required <- c("Seurat", "SeuratObject", "Matrix", "Biobase", "monocle", "VGAM", "future", "igraph", "ggplot2", "dplyr", 
        "tidyr", "tibble", "scatterpie", "scales", "mgcv")
    missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
    if (length(missing) > 0L) {
        stop("Missing package(s): ", paste(missing, collapse = ", "), ". This script does not install packages automatically.")
    }
    suppressPackageStartupMessages({
        library(Seurat)
        library(SeuratObject)
        library(Matrix)
        library(Biobase)
        library(monocle)
        library(VGAM)
        library(future)
        library(igraph)
        library(ggplot2)
        library(dplyr)
        library(tidyr)
        library(tibble)
        library(scatterpie)
        library(scales)
        library(mgcv)
    })
    future::plan(future::sequential)
    monocle_version <- as.character(utils::packageVersion("monocle"))
    monocle_path <- normalizePath(find.package("monocle"), winslash = "/")
    message("Monocle version: ", monocle_version)
    message("Monocle installed path: ", monocle_path)
    message("R libraries: ", paste(.libPaths(), collapse = " | "))
    if (!identical(monocle_version, EXPECTED_MONOCLE_VERSION)) {
        stop("Expected modified monocle ", EXPECTED_MONOCLE_VERSION, ", but found ", monocle_version, ".")
    }
    write_csv <- function(x, file) {
        utils::write.csv(x, file, row.names = FALSE, fileEncoding = "UTF-8")
    }
    save_plot <- function(p, stem, width = 7.2, height = 5.8) {
        ggplot2::ggsave(file.path(FIGURE_DIR, paste0(stem, ".pdf")), p, width = width, height = height, units = "in", device = "pdf", 
            bg = "white")
        ggplot2::ggsave(file.path(FIGURE_DIR, paste0(stem, ".png")), p, width = width, height = height, units = "in", dpi = 320, 
            bg = "white")
    }
    template_theme <- function(base_size = 14) {
        ggplot2::theme_bw(base_size = base_size) + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", 
            size = 15), panel.grid = ggplot2::element_blank(), axis.title = ggplot2::element_text(face = "bold", size = 13), 
            axis.text = ggplot2::element_text(color = "black", size = 11), legend.title = ggplot2::element_text(face = "bold", 
                size = 11), legend.text = ggplot2::element_text(size = 10))
    }
    trajectory_plot <- function(cds, color_by, title, discrete_colors = NULL, continuous_colors = NULL, cell_size = 0.65) {
        p <- monocle::plot_cell_trajectory(cds, color_by = color_by, cell_size = cell_size, show_branch_points = TRUE) + 
            ggplot2::ggtitle(title) + template_theme()
        if (!is.null(discrete_colors)) {
            p <- p + ggplot2::scale_color_manual(values = discrete_colors, drop = FALSE)
        }
        if (!is.null(continuous_colors)) {
            p <- p + ggplot2::scale_color_gradientn(colours = continuous_colors, name = color_by)
        }
        p
    }
    map_gene <- function(features, gene) {
        idx <- which(toupper(features) == toupper(gene))
        if (length(idx) != 1L) 
            stop("Cannot uniquely match gene: ", gene)
        features[idx]
    }
    get_sct_data <- function(object) {
        if (!"SCT" %in% SeuratObject::Assays(object)) 
            stop("Final object has no SCT assay.")
        if (!"data" %in% SeuratObject::Layers(object[["SCT"]])) {
            stop("Final object SCT assay has no data layer.")
        }
        SeuratObject::LayerData(object[["SCT"]], layer = "data")
    }
    get_count_layers <- function(object, assay_name, target_cells) {
        if (!assay_name %in% SeuratObject::Assays(object)) {
            stop("Assay not found: ", assay_name)
        }
        assay <- object[[assay_name]]
        count_layers <- SeuratObject::Layers(assay)
        count_layers <- count_layers[grepl("^counts($|\\.)", count_layers)]
        if (length(count_layers) == 0L) 
            stop("No counts layer in assay: ", assay_name)
        mats <- lapply(count_layers, function(layer_name) {
            x <- SeuratObject::LayerData(assay, layer = layer_name)
            keep <- intersect(target_cells, colnames(x))
            if (length(keep) == 0L) 
                return(NULL)
            x[, keep, drop = FALSE]
        })
        mats <- Filter(Negate(is.null), mats)
        if (length(mats) == 0L) 
            stop("No target cells found in count layers.")
        common_genes <- Reduce(intersect, lapply(mats, rownames))
        if (length(common_genes) == 0L) 
            stop("Count layers have no common genes.")
        mats <- lapply(mats, function(x) x[common_genes, , drop = FALSE])
        counts <- do.call(cbind, mats)
        if (anyDuplicated(colnames(counts))) 
            stop("Duplicated barcodes in counts.")
        missing_cells <- setdiff(target_cells, colnames(counts))
        if (length(missing_cells) > 0L) {
            stop(length(missing_cells), " final fibroblasts are missing from counts.")
        }
        counts <- methods::as(counts[, target_cells, drop = FALSE], "dgCMatrix")
        if (length(counts@x) && any(!is.finite(counts@x) | counts@x < 0)) {
            stop("Counts contain invalid values.")
        }
        if (length(counts@x) && any(abs(counts@x - round(counts@x)) > 1e-08)) {
            stop("Selected counts are not integer-like; normalized data will not be used.")
        }
        counts
    }
    get_sample_field <- function(meta) {
        hit <- intersect(c("sample_id", "orig.ident", "sample", "donor_id"), colnames(meta))
        if (length(hit)) 
            hit[1]
        else NA
    }
    order_cells_on_principal_graph <- function(cds, root_vertex = NULL, cell_state_override = NULL) {
        mst <- monocle::minSpanningTree(cds)
        Y <- monocle::reducedDimK(cds)
        Z <- monocle::reducedDimS(cds)
        if (is.null(mst) || ncol(Y) < 2L || ncol(Z) != ncol(cds)) {
            stop("DDRTree principal graph is unavailable or malformed.")
        }
        closest <- cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_closest_vertex
        if (is.null(closest)) {
            cds <- monocle:::findNearestPointOnMST(cds)
            closest <- cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_closest_vertex
        }
        cell_names <- colnames(cds)
        vertex_names <- colnames(Y)
        closest_idx <- as.integer(closest[cell_names, 1])
        if (length(closest_idx) != length(cell_names) || any(!is.finite(closest_idx)) || any(closest_idx < 1L | closest_idx > 
            length(vertex_names))) {
            stop("Cell-to-principal-vertex mapping is invalid.")
        }
        if (is.null(root_vertex)) {
            diameter_path <- igraph::get_diameter(mst, directed = FALSE, weights = igraph::E(mst)$weight)
            root_vertex <- igraph::as_ids(diameter_path)[1]
        }
        root_vertex <- as.character(root_vertex)[1]
        if (!root_vertex %in% vertex_names) 
            stop("Root vertex is absent from DDRTree.")
        vertex_order <- monocle:::extract_ddrtree_ordering(cds, root_vertex)
        vertex_pt <- stats::setNames(as.numeric(vertex_order$pseudo_time), rownames(vertex_order))[vertex_names]
        vertex_state <- stats::setNames(as.character(vertex_order$cell_state), rownames(vertex_order))[vertex_names]
        if (any(!is.finite(vertex_pt)) || anyNA(vertex_state)) {
            stop("Principal-vertex ordering returned invalid values.")
        }
        cell_pt <- numeric(length(cell_names))
        for (v_idx in sort(unique(closest_idx))) {
            idx <- which(closest_idx == v_idx)
            v_name <- vertex_names[v_idx]
            neighbors <- igraph::as_ids(igraph::neighbors(mst, v_name, mode = "all"))
            if (length(neighbors) == 0L) {
                cell_pt[idx] <- vertex_pt[v_name]
                next
            }
            z <- Z[, idx, drop = FALSE]
            y_v <- Y[, v_name]
            best_d2 <- rep(Inf, length(idx))
            best_pt <- rep(vertex_pt[v_name], length(idx))
            for (u_name in neighbors) {
                edge_vector <- Y[, u_name] - y_v
                edge_length_sq <- sum(edge_vector^2)
                if (!is.finite(edge_length_sq) || edge_length_sq <= 0) 
                  next
                centered <- sweep(z, 1L, y_v, FUN = "-")
                edge_fraction <- as.numeric(crossprod(edge_vector, centered))/edge_length_sq
                edge_fraction <- pmin(1, pmax(0, edge_fraction))
                projected <- y_v + tcrossprod(edge_vector, edge_fraction)
                distance_sq <- colSums((z - projected)^2)
                projected_pt <- (1 - edge_fraction) * vertex_pt[v_name] + edge_fraction * vertex_pt[u_name]
                better <- distance_sq < best_d2
                best_d2[better] <- distance_sq[better]
                best_pt[better] <- projected_pt[better]
            }
            cell_pt[idx] <- best_pt
        }
        if (any(!is.finite(cell_pt)) || length(cell_pt) != ncol(cds)) {
            stop("Memory-safe cell pseudotime contains invalid values.")
        }
        cell_pt <- cell_pt - min(cell_pt)
        names(cell_pt) <- cell_names
        if (is.null(cell_state_override)) {
            cell_state <- vertex_state[vertex_names[closest_idx]]
        }
        else {
            if (length(cell_state_override) != ncol(cds)) {
                stop("State override length does not match the cells.")
            }
            cell_state <- as.character(cell_state_override)
        }
        Biobase::pData(cds)$Pseudotime <- unname(cell_pt[cell_names])
        Biobase::pData(cds)$State <- factor(cell_state)
        cds@auxOrderingData[["DDRTree"]]$root_cell <- root_vertex
        cds@auxOrderingData[["DDRTree"]]$branch_points <- igraph::as_ids(igraph::V(mst)[igraph::degree(mst) > 2])
        cds@auxOrderingData[["DDRTree"]]$principal_vertex_pseudotime <- vertex_pt
        cds@auxOrderingData[["DDRTree"]]$principal_vertex_state <- vertex_state
        cds@auxOrderingData[["DDRTree"]]$ordering_method <- "principal_graph_edge_projection_memory_safe"
        cds
    }
    choose_root_vertex_for_state <- function(cds, root_state) {
        aux <- cds@auxOrderingData[["DDRTree"]]
        vertex_state <- aux$principal_vertex_state
        vertex_pt <- aux$principal_vertex_pseudotime
        current_root <- as.character(aux$root_cell)[1]
        candidates <- names(vertex_state)[vertex_state == as.character(root_state)]
        if (length(candidates) == 0L) {
            stop("No DDRTree principal vertex belongs to selected root State: ", root_state)
        }
        if (current_root %in% candidates) {
            candidates[which.min(vertex_pt[candidates])]
        }
        else {
            candidates[which.max(vertex_pt[candidates])]
        }
    }
    for (f in c(FINAL_FIBRO_RDS, WHOLE_ATLAS_RDS, OXPHOS_SCORE_CSV)) {
        if (!file.exists(f)) 
            stop("Input file not found: ", f)
    }
    message("Reading locked final fibroblast object ...")
    fib <- readRDS(FINAL_FIBRO_RDS)
    message("Reading completed whole atlas ...")
    atlas <- readRDS(WHOLE_ATLAS_RDS)
    if (!inherits(fib, "Seurat") || !inherits(atlas, "Seurat")) {
        stop("Both RDS files must contain Seurat objects.")
    }
    if (ncol(fib) != EXPECTED_FIBRO_CELLS) {
        stop("Final fibroblast object has ", ncol(fib), " cells; expected ", EXPECTED_FIBRO_CELLS, ".")
    }
    if (ncol(atlas) != EXPECTED_ATLAS_CELLS) {
        stop("Whole atlas has ", ncol(atlas), " cells; expected ", EXPECTED_ATLAS_CELLS, ".")
    }
    if (!all(c(GROUP_FIELD, SUBTYPE_FIELD) %in% colnames(fib@meta.data))) {
        stop("Final object lacks Group or fibroblast_subtype_published metadata.")
    }
    if (!GROUP_FIELD %in% colnames(atlas@meta.data)) 
        stop("Whole atlas lacks Group.")
    if (anyDuplicated(colnames(fib)) || anyDuplicated(colnames(atlas))) {
        stop("Duplicated cell barcodes detected.")
    }
    fib_cells <- colnames(fib)
    if (!all(fib_cells %in% colnames(atlas))) {
        stop(sum(!fib_cells %in% colnames(atlas)), " final fibroblast barcodes are absent from the whole atlas.")
    }
    observed_subtypes <- table(as.character(fib@meta.data[[SUBTYPE_FIELD]]))
    if (!setequal(names(observed_subtypes), SUBTYPE_LEVELS) || !identical(as.integer(observed_subtypes[SUBTYPE_LEVELS]), 
        as.integer(EXPECTED_SUBTYPE_COUNTS[SUBTYPE_LEVELS]))) {
        stop("Locked seven-subtype labels or counts have changed.")
    }
    fib_group <- as.character(fib@meta.data[fib_cells, GROUP_FIELD])
    atlas_group <- as.character(atlas@meta.data[fib_cells, GROUP_FIELD])
    if (!identical(fib_group, atlas_group)) 
        stop("Group labels differ between RDS files.")
    if (!setequal(unique(fib_group), CONDITION_LEVELS)) {
        stop("Group must contain exactly Control and OA.")
    }
    sct_data <- get_sct_data(fib)
    ndufs3_sct_gene <- map_gene(rownames(sct_data), TARGET_GENE)
    ndufs3_sct <- as.numeric(sct_data[ndufs3_sct_gene, fib_cells])
    names(ndufs3_sct) <- fib_cells
    if (any(!is.finite(ndufs3_sct)) || any(ndufs3_sct < 0)) {
        stop("NDUFS3 SCT/data contains invalid values.")
    }
    oxphos <- utils::read.csv(OXPHOS_SCORE_CSV, check.names = FALSE)
    needed_score_cols <- c("cell_id", "score_consensus_z")
    if (!all(needed_score_cols %in% colnames(oxphos))) {
        stop("OXPHOS score CSV lacks: ", paste(setdiff(needed_score_cols, colnames(oxphos)), collapse = ", "))
    }
    if (anyDuplicated(oxphos$cell_id)) 
        stop("OXPHOS score CSV has duplicated cell_id.")
    if (!setequal(oxphos$cell_id, fib_cells)) {
        stop("OXPHOS score cell IDs do not exactly match the 31,004 final fibroblasts.")
    }
    oxphos <- oxphos[match(fib_cells, oxphos$cell_id), , drop = FALSE]
    if (any(!is.finite(oxphos$score_consensus_z))) {
        stop("OXPHOS score_consensus_z contains non-finite values.")
    }
    atlas_has_rna_counts <- FALSE
    if ("RNA" %in% SeuratObject::Assays(atlas)) {
        atlas_has_rna_counts <- any(grepl("^counts($|\\.)", SeuratObject::Layers(atlas[["RNA"]])))
    }
    if (atlas_has_rna_counts) {
        message("Using whole-atlas RNA raw counts.")
        raw_counts <- get_count_layers(atlas, "RNA", fib_cells)
        count_source <- "whole-atlas RNA raw counts"
    }
    else {
        message("Whole atlas has no RNA counts; using final-fibroblast SCT corrected counts.")
        raw_counts <- get_count_layers(fib, "SCT", fib_cells)
        count_source <- "final-fibroblast SCT corrected counts"
    }
    ndufs3_count_gene <- map_gene(rownames(raw_counts), TARGET_GENE)
    sample_field <- get_sample_field(fib@meta.data)
    sample_values <- if (is.na(sample_field)) {
        rep("sample_not_available", length(fib_cells))
    }
    else {
        x <- as.character(fib@meta.data[fib_cells, sample_field])
        x[is.na(x) | x == ""] <- "sample_not_available"
        x
    }
    pd <- fib@meta.data[fib_cells, , drop = FALSE]
    pd$Group <- factor(fib_group, levels = CONDITION_LEVELS)
    pd$fibroblast_subtype_published <- factor(as.character(pd[[SUBTYPE_FIELD]]), levels = SUBTYPE_LEVELS)
    pd$sample_for_plot <- factor(sample_values, levels = unique(sample_values))
    pd$NDUFS3_status <- factor(ifelse(ndufs3_sct > 0, "NDUFS3+", "NDUFS3-"), levels = c("NDUFS3-", "NDUFS3+"))
    pd$OXPHOS_consensus_z <- oxphos$score_consensus_z
    pd$PRG4_NDUFS3_state <- "Other fibroblasts"
    is_prg4 <- as.character(pd$fibroblast_subtype_published) == PRG4_SUBTYPE
    pd$PRG4_NDUFS3_state[is_prg4 & pd$NDUFS3_status == "NDUFS3+"] <- "PRG4+ / NDUFS3+"
    pd$PRG4_NDUFS3_state[is_prg4 & pd$NDUFS3_status == "NDUFS3-"] <- "PRG4+ / NDUFS3-"
    pd$PRG4_NDUFS3_state <- factor(pd$PRG4_NDUFS3_state, levels = c("Other fibroblasts", "PRG4+ / NDUFS3-", "PRG4+ / NDUFS3+"))
    if (!identical(colnames(raw_counts), rownames(pd))) 
        stop("Counts and metadata misaligned.")
    qc <- data.frame(metric = c("whole_atlas_cells", "locked_final_fibroblast_cells", "exact_barcode_matches", "unmatched_cells", 
        "count_source", "count_genes", "NDUFS3_positive", "NDUFS3_negative", "OXPHOS_score_matches", "monocle_version", "monocle_path"), 
        value = c(ncol(atlas), ncol(fib), length(fib_cells), 0, count_source, nrow(raw_counts), sum(pd$NDUFS3_status == "NDUFS3+"), 
            sum(pd$NDUFS3_status == "NDUFS3-"), nrow(oxphos), monocle_version, monocle_path))
    write_csv(qc, file.path(TABLE_DIR, "01_input_and_package_QC.csv"))
    write_csv(as.data.frame(table(Group = pd$Group, fibroblast_subtype = pd$fibroblast_subtype_published, NDUFS3_status = pd$NDUFS3_status)), 
        file.path(TABLE_DIR, "02_cell_counts_Group_subtype_NDUFS3.csv"))
    print(qc)
    rm(atlas, fib, sct_data, oxphos)
    invisible(gc(full = TRUE))
    marker_file <- file.path(TABLE_DIR, "03_locked_subtype_markers.csv")
    ordering_file <- file.path(TABLE_DIR, "04_ordering_genes.csv")
    if (REUSE_MARKER_TABLES && file.exists(marker_file) && file.exists(ordering_file)) {
        message("Reusing completed marker and ordering-gene tables.")
        subtype_markers <- utils::read.csv(marker_file, check.names = FALSE)
        ordering_table <- utils::read.csv(ordering_file, check.names = FALSE)
        ordering_genes <- unique(ordering_table$gene)
        ordering_genes <- ordering_genes[ordering_genes %in% rownames(raw_counts)]
        ordering_genes <- head(ordering_genes, MAX_ORDERING_GENES)
    }
    else {
        message("Finding markers among the seven locked fibroblast subtypes ...")
        marker_object <- Seurat::CreateSeuratObject(counts = raw_counts, meta.data = pd, assay = "RNA", project = "GSE216651_locked_fibroblasts")
        marker_object <- Seurat::NormalizeData(marker_object, verbose = FALSE)
        Seurat::Idents(marker_object) <- "fibroblast_subtype_published"
        subtype_markers <- Seurat::FindAllMarkers(marker_object, assay = "RNA", only.pos = FALSE, min.pct = MIN_PCT, logfc.threshold = LOGFC_THRESHOLD, 
            test.use = "wilcox", verbose = TRUE)
        if (!"gene" %in% colnames(subtype_markers)) {
            subtype_markers <- tibble::rownames_to_column(subtype_markers, "gene")
        }
        fc_col <- intersect(c("avg_log2FC", "avg_logFC"), colnames(subtype_markers))[1]
        if (is.na(fc_col)) 
            stop("FindAllMarkers result has no logFC column.")
        write_csv(subtype_markers, marker_file)
        ordering_table <- subtype_markers %>% dplyr::filter(is.finite(p_val_adj), p_val_adj < ORDERING_FDR) %>% dplyr::mutate(abs_logFC = abs(.data[[fc_col]])) %>% 
            dplyr::group_by(gene) %>% dplyr::summarise(best_FDR = min(p_val_adj), max_abs_logFC = max(abs_logFC), .groups = "drop") %>% 
            dplyr::arrange(best_FDR, dplyr::desc(max_abs_logFC)) %>% dplyr::filter(gene %in% rownames(raw_counts)) %>% dplyr::slice_head(n = MAX_ORDERING_GENES) %>% 
            dplyr::mutate(ordering_rank = dplyr::row_number())
        ordering_genes <- ordering_table$gene
        write_csv(ordering_table, ordering_file)
        rm(marker_object)
        invisible(gc(full = TRUE))
    }
    if (length(ordering_genes) < MIN_ORDERING_GENES) {
        stop("Only ", length(ordering_genes), " usable ordering genes were retained.")
    }
    message("Ordering genes: ", length(ordering_genes), "; NDUFS3 included: ", TARGET_GENE %in% ordering_genes)
    pre_ddrtree_file <- file.path(RDS_DIR, "01_Fibroblasts_before_DDRTree.rds")
    if (REUSE_SAVED_CDS && file.exists(pre_ddrtree_file)) {
        message("Reusing saved CellDataSet before DDRTree: ", pre_ddrtree_file)
        cds <- readRDS(pre_ddrtree_file)
        if (!inherits(cds, "CellDataSet") || ncol(cds) != EXPECTED_FIBRO_CELLS || !identical(colnames(cds), rownames(pd))) {
            stop("Saved pre-DDRTree CellDataSet does not match the current 31,004 cells.")
        }
    }
    else {
        fd <- data.frame(gene_short_name = rownames(raw_counts), row.names = rownames(raw_counts), stringsAsFactors = FALSE)
        cds <- monocle::newCellDataSet(raw_counts, phenoData = new("AnnotatedDataFrame", data = pd), featureData = new("AnnotatedDataFrame", 
            data = fd), lowerDetectionLimit = 0.5, expressionFamily = VGAM::negbinomial.size())
        cds <- BiocGenerics::estimateSizeFactors(cds)
        cds <- BiocGenerics::estimateDispersions(cds)
        cds <- monocle::detectGenes(cds, min_expr = 0.1)
        cds <- monocle::setOrderingFilter(cds, ordering_genes)
        saveRDS(cds, pre_ddrtree_file)
    }
    pdf(file.path(FIGURE_DIR, "00_ordering_gene_dispersion.pdf"), 7, 6)
    tryCatch(print(monocle::plot_ordering_genes(cds)), finally = dev.off())
    png(file.path(FIGURE_DIR, "00_ordering_gene_dispersion.png"), width = 7, height = 6, units = "in", res = 320, bg = "white")
    tryCatch(print(monocle::plot_ordering_genes(cds)), finally = dev.off())
    rm(raw_counts)
    if (exists("fd")) 
        rm(fd)
    invisible(gc(full = TRUE))
    ddrtree_file <- file.path(RDS_DIR, "01B_Fibroblasts_DDRTree_before_ordering.rds")
    if (REUSE_SAVED_CDS && file.exists(ddrtree_file)) {
        message("Reusing completed DDRTree principal graph: ", ddrtree_file)
        cds <- readRDS(ddrtree_file)
        if (!inherits(cds, "CellDataSet") || ncol(cds) != EXPECTED_FIBRO_CELLS || is.null(monocle::minSpanningTree(cds))) {
            stop("Saved DDRTree CellDataSet is invalid.")
        }
    }
    else {
        message("Running DDRTree on all 31,004 locked final fibroblasts ...")
        cds <- monocle::reduceDimension(cds, max_components = 2, method = "DDRTree")
        saveRDS(cds, ddrtree_file)
    }
    message("Ordering all cells by memory-safe DDRTree principal-graph projection ...")
    cds <- order_cells_on_principal_graph(cds)
    initial_meta <- Biobase::pData(cds)
    state_group <- as.matrix(table(initial_meta$State, initial_meta$Group))
    if (!all(CONDITION_LEVELS %in% colnames(state_group))) {
        stop("Both Control and OA are required for root selection.")
    }
    root_table <- data.frame(State = rownames(state_group), Control_n = as.numeric(state_group[, "Control"]), OA_n = as.numeric(state_group[, 
        "OA"]))
    root_table$total_n <- root_table$Control_n + root_table$OA_n
    root_table$Control_proportion <- root_table$Control_n/root_table$total_n
    root_table$OA_proportion <- root_table$OA_n/root_table$total_n
    root_table <- root_table[order(-root_table$Control_proportion, -root_table$Control_n, suppressWarnings(as.numeric(root_table$State))), 
        ]
    root_table$selected_root <- FALSE
    root_table$selected_root[1] <- TRUE
    root_state <- suppressWarnings(as.numeric(root_table$State[1]))
    if (!is.finite(root_state)) 
        stop("Selected root State is not numeric.")
    write_csv(root_table, file.path(TABLE_DIR, "05_root_State_Control_proportion.csv"))
    initial_cell_state <- Biobase::pData(cds)$State
    root_vertex <- choose_root_vertex_for_state(cds, root_state)
    message("Selected root principal vertex: ", root_vertex)
    cds <- order_cells_on_principal_graph(cds, root_vertex = root_vertex, cell_state_override = initial_cell_state)
    size_factors <- BiocGenerics::sizeFactors(cds)
    if (any(!is.finite(size_factors) | size_factors <= 0)) 
        stop("Invalid size factors.")
    ndufs3_raw <- as.numeric(Biobase::exprs(cds)[ndufs3_count_gene, colnames(cds)])
    Biobase::pData(cds)$NDUFS3_expression_log1p_normalized <- log1p(ndufs3_raw/size_factors)
    saveRDS(cds, file.path(RDS_DIR, "02_Fibroblasts_ordered.rds"))
    meta <- Biobase::pData(cds) %>% as.data.frame() %>% tibble::rownames_to_column("cell")
    write_csv(meta, file.path(TABLE_DIR, "06_cell_level_State_Pseudotime_metadata.csv"))
    write_state_table <- function(field, filename) {
        out <- meta %>% dplyr::count(State, .data[[field]], name = "n_cells") %>% dplyr::group_by(State) %>% dplyr::mutate(proportion_within_State = n_cells/sum(n_cells)) %>% 
            dplyr::ungroup()
        write_csv(out, file.path(TABLE_DIR, filename))
    }
    write_state_table("Group", "07_State_composition_by_Group.csv")
    write_state_table("fibroblast_subtype_published", "08_State_composition_by_subtype.csv")
    write_state_table("NDUFS3_status", "09_State_composition_by_NDUFS3.csv")
    plot_state_composition <- function(field, colors, title, legend_title, stem) {
        plot_data <- meta %>% dplyr::count(State, .data[[field]], name = "n_cells") %>% dplyr::group_by(State) %>% dplyr::mutate(proportion = n_cells/sum(n_cells)) %>% 
            dplyr::ungroup()
        p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = factor(State), y = proportion, fill = .data[[field]])) + ggplot2::geom_col(width = 0.78, 
            color = "white", linewidth = 0.25) + ggplot2::scale_fill_manual(values = colors, drop = FALSE) + ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1)) + 
            ggplot2::labs(x = "State", y = "Proportion within State", fill = legend_title, title = title) + template_theme()
        save_plot(p, stem, 7.4, 5.6)
    }
    p_subtype <- trajectory_plot(cds, "fibroblast_subtype_published", "Fibroblast trajectory by locked subtype", discrete_colors = SUBTYPE_COLORS)
    save_plot(p_subtype, "01_trajectory_locked_7_subtypes", 8.2, 6.3)
    p_state <- trajectory_plot(cds, "State", "Fibroblast trajectory by State")
    save_plot(p_state, "02_trajectory_State")
    p_pseudotime <- trajectory_plot(cds, "Pseudotime", "Fibroblast trajectory by Pseudotime", continuous_colors = PSEUDOTIME_COLORS)
    save_plot(p_pseudotime, "03_trajectory_Pseudotime")
    p_group <- trajectory_plot(cds, "Group", "Fibroblast trajectory by Group", discrete_colors = GROUP_COLORS)
    save_plot(p_group, "04_trajectory_Control_OA")
    p_status <- trajectory_plot(cds, "NDUFS3_status", "NDUFS3 transcript states on the fibroblast trajectory", discrete_colors = NDUFS3_COLORS)
    save_plot(p_status, "05_trajectory_NDUFS3_positive_negative")
    p_expression <- trajectory_plot(cds, "NDUFS3_expression_log1p_normalized", "NDUFS3 expression on the fibroblast trajectory", 
        continuous_colors = c("#F7F7F7", "#FDB863", "#B2182B"))
    save_plot(p_expression, "06_trajectory_NDUFS3_expression")
    prg4_colors <- c(`Other fibroblasts` = "#D9D9D9", `PRG4+ / NDUFS3-` = NDUFS3_COLORS[["NDUFS3-"]], `PRG4+ / NDUFS3+` = NDUFS3_COLORS[["NDUFS3+"]])
    p_prg4 <- trajectory_plot(cds, "PRG4_NDUFS3_state", "PRG4+ NDUFS3 states on the fibroblast trajectory", discrete_colors = prg4_colors)
    save_plot(p_prg4, "07_trajectory_PRG4_NDUFS3_states")
    p_oxphos <- trajectory_plot(cds, "OXPHOS_consensus_z", "Oxidative phosphorylation score on the trajectory", continuous_colors = OXPHOS_COLORS)
    save_plot(p_oxphos, "08_trajectory_OXPHOS_consensus_score")
    trend_data <- meta %>% dplyr::select(cell, Pseudotime, State, Group, fibroblast_subtype_published, NDUFS3_status, NDUFS3_expression_log1p_normalized, 
        OXPHOS_consensus_z)
    write_csv(trend_data, file.path(TABLE_DIR, "10_NDUFS3_OXPHOS_along_Pseudotime.csv"))
    p_ndufs3_trend <- ggplot2::ggplot(trend_data, ggplot2::aes(Pseudotime, NDUFS3_expression_log1p_normalized, color = Group)) + 
        ggplot2::geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), se = TRUE, linewidth = 1, alpha = 0.16) + ggplot2::scale_color_manual(values = GROUP_COLORS, 
        drop = FALSE) + ggplot2::labs(x = "Pseudotime", y = "NDUFS3 expression\n(log1p normalized counts)", color = NULL, 
        title = "NDUFS3 expression along fibroblast pseudotime") + template_theme()
    save_plot(p_ndufs3_trend, "09_NDUFS3_expression_along_Pseudotime", 7.2, 5.5)
    p_oxphos_trend <- ggplot2::ggplot(trend_data, ggplot2::aes(Pseudotime, OXPHOS_consensus_z, color = Group)) + ggplot2::geom_smooth(method = "gam", 
        formula = y ~ s(x, bs = "cs"), se = TRUE, linewidth = 1, alpha = 0.16) + ggplot2::scale_color_manual(values = GROUP_COLORS, 
        drop = FALSE) + ggplot2::labs(x = "Pseudotime", y = "OXPHOS five-method consensus score (z)", color = NULL, title = "Oxidative phosphorylation score along pseudotime") + 
        template_theme()
    save_plot(p_oxphos_trend, "10_OXPHOS_score_along_Pseudotime", 7.2, 5.5)
    sample_colors <- setNames((scales::hue_pal(l = 65, c = 90))(nlevels(Biobase::pData(cds)$sample_for_plot)), levels(Biobase::pData(cds)$sample_for_plot))
    p_sample <- trajectory_plot(cds, "sample_for_plot", "Fibroblast trajectory by sample (QC)", discrete_colors = sample_colors, 
        cell_size = 0.55)
    save_plot(p_sample, "11_trajectory_sample_QC", 8.2, 6.3)
    plot_state_composition("Group", GROUP_COLORS, "Control and OA composition within each State", NULL, "15_State_composition_Control_OA")
    plot_state_composition("fibroblast_subtype_published", SUBTYPE_COLORS, "Locked fibroblast-subtype composition within each State", 
        "Subtype", "16_State_composition_locked_7_subtypes")
    plot_state_composition("NDUFS3_status", NDUFS3_COLORS, "NDUFS3 transcript-state composition within each State", NULL, 
        "17_State_composition_NDUFS3_positive_negative")
    rd <- monocle::reducedDimS(cds)
    coord <- as.data.frame(t(rd[1:2, , drop = FALSE]))
    colnames(coord) <- c("Comp1", "Comp2")
    coord$cell <- rownames(coord)
    trajectory_df <- dplyr::left_join(coord, meta, by = "cell")
    make_pie_data <- function(df) {
        state_prop <- df %>% dplyr::count(State, NDUFS3_status, name = "n") %>% tidyr::pivot_wider(names_from = NDUFS3_status, 
            values_from = n, values_fill = 0)
        if (!"NDUFS3+" %in% colnames(state_prop)) 
            state_prop[["NDUFS3+"]] <- 0
        if (!"NDUFS3-" %in% colnames(state_prop)) 
            state_prop[["NDUFS3-"]] <- 0
        state_prop <- state_prop %>% dplyr::mutate(total = .data[["NDUFS3+"]] + .data[["NDUFS3-"]], positive_ratio = .data[["NDUFS3+"]]/total)
        centers <- df %>% dplyr::group_by(State) %>% dplyr::summarise(Comp1 = median(Comp1, na.rm = TRUE), Comp2 = median(Comp2, 
            na.rm = TRUE), .groups = "drop")
        pie <- dplyr::left_join(centers, state_prop, by = "State")
        pie$r <- min(diff(range(trajectory_df$Comp1, na.rm = TRUE)), diff(range(trajectory_df$Comp2, na.rm = TRUE))) * 0.06
        pie
    }
    add_state_pies <- function(base_plot, pie_data, title) {
        base_plot + scatterpie::geom_scatterpie(data = pie_data, ggplot2::aes(x = Comp1, y = Comp2, r = r), cols = c("NDUFS3+", 
            "NDUFS3-"), inherit.aes = FALSE, color = "black", alpha = 0.95, linewidth = 0.5) + ggplot2::geom_text(data = pie_data, 
            ggplot2::aes(x = Comp1, y = Comp2, label = State), inherit.aes = FALSE, size = 4, fontface = "bold", color = "black") + 
            ggplot2::scale_fill_manual(values = NDUFS3_COLORS, name = "group") + ggplot2::ggtitle(title) + ggplot2::coord_cartesian(clip = "off") + 
            template_theme()
    }
    p_pie_base <- monocle::plot_cell_trajectory(cds, color_by = "Pseudotime", cell_size = 1.2) + ggplot2::scale_color_gradientn(colours = PSEUDOTIME_COLORS, 
        name = "Pseudotime") + template_theme()
    pie_all <- make_pie_data(trajectory_df)
    write_csv(pie_all, file.path(TABLE_DIR, "11_State_pie_NDUFS3_all_fibroblasts.csv"))
    p_pie_all <- add_state_pies(p_pie_base, pie_all, "NDUFS3 positive proportion in each State")
    save_plot(p_pie_all, "12_State_pie_NDUFS3_all_fibroblasts", 8.5, 6.5)
    trajectory_prg4 <- trajectory_df %>% dplyr::filter(as.character(fibroblast_subtype_published) == PRG4_SUBTYPE)
    pie_prg4 <- make_pie_data(trajectory_prg4)
    write_csv(pie_prg4, file.path(TABLE_DIR, "12_State_pie_NDUFS3_PRG4_only.csv"))
    p_pie_prg4 <- add_state_pies(p_pie_base, pie_prg4, "NDUFS3 positive proportion among PRG4+ cells in each State")
    save_plot(p_pie_prg4, "13_State_pie_NDUFS3_PRG4_only", 8.5, 6.5)
    invisible(NULL)
    write_csv(parameters, file.path(TABLE_DIR, "00_analysis_parameters.csv"))
    saveRDS(list(cds = cds, ordering_genes = ordering_genes, subtype_markers = subtype_markers, root_table = root_table), 
        file.path(RDS_DIR, "03_Fibroblasts_Monocle2_final_results.rds"))
    capture.output(sessionInfo(), file = file.path(LOG_DIR, "sessionInfo.txt"))
    message("Analysis completed: ", format(Sys.time()))
    message("Results: ", OUTPUT_DIR)
})

