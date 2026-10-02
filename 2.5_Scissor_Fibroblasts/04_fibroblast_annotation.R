args <- commandArgs(trailingOnly = FALSE)

script_file <- sub("^--file=", "", grep("^--file=", args, value = TRUE)[1])

code_root <- dirname(dirname(normalizePath(script_file, mustWork = TRUE)))

source(file.path(code_root, "00_configuration.R"))

gc()

suppressPackageStartupMessages({
    library(Seurat)
    library(ggplot2)
})

options(stringsAsFactors = FALSE)

set.seed(20260731)

input_rds <- paste0(paste0(file.path(input_root, "singlecell"), "/"), "20_Fibroblast_secondary_clustering/", "03_clean_recluster_after_cluster8_removal_corrected/", 
    "01_rds/", "GSE216651_Fibroblasts_without_cluster8_PC20_res0.5.rds")

output_dir <- paste0(paste0(file.path(input_root, "singlecell"), "/"), "20_Fibroblast_secondary_clustering/", "07_final_fibroblast_annotation")

rds_dir <- file.path(output_dir, "01_rds")

table_dir <- file.path(output_dir, "02_tables")

figure_dir <- file.path(output_dir, "03_figures")

log_dir <- file.path(output_dir, "04_logs")

for (x in c(rds_dir, table_dir, figure_dir, log_dir)) {
    dir.create(x, recursive = TRUE, showWarnings = FALSE)
}

log_file <- file.path(log_dir, "20H_final_annotation_run.log")

log_message <- function(...) {
    txt <- paste0(...)
    message(txt)
    cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", txt, "\n", file = log_file, append = TRUE, sep = "")
}

invisible(NULL)

if (!file.exists(input_rds)) {
    stop("Input RDS does not exist: ", input_rds)
}

fib <- readRDS(input_rds)

cluster_column <- "seurat_clusters"

sample_column <- "sample_id"

if (!"SCT" %in% Assays(fib)) {
    stop("The SCT assay required for the marker DotPlot is absent. ", "Available assays: ", paste(Assays(fib), collapse = ", "))
}

DefaultAssay(fib) <- "SCT"

required_metadata <- c(cluster_column, sample_column)

missing_metadata <- setdiff(required_metadata, colnames(fib[[]]))

if (length(missing_metadata) > 0) {
    stop("Required metadata column(s) missing: ", paste(missing_metadata, collapse = ", "))
}

if (ncol(fib) == 0) {
    stop("The input Seurat object contains no cells.")
}

original_cluster <- as.character(fib[[cluster_column]][, 1])

if (anyNA(original_cluster) || any(original_cluster == "")) {
    stop("The original cluster column contains missing or empty labels.")
}

expected_original_clusters <- as.character(0:9)

observed_original_clusters <- sort(unique(original_cluster))

unexpected_clusters <- setdiff(observed_original_clusters, expected_original_clusters)

missing_clusters <- setdiff(expected_original_clusters, observed_original_clusters)

if (length(unexpected_clusters) > 0) {
    stop("Unexpected original cluster label(s): ", paste(unexpected_clusters, collapse = ", "))
}

if (length(missing_clusters) > 0) {
    stop("Expected original cluster(s) absent: ", paste(missing_clusters, collapse = ", "), ". Confirm that the correct RDS was selected.")
}

sample_values <- as.character(fib[[sample_column]][, 1])

if (anyNA(sample_values) || any(sample_values == "")) {
    stop("The sample_id column contains missing or empty values.")
}

log_message("Cells loaded: ", ncol(fib))

log_message("Original clusters: ", paste(observed_original_clusters, collapse = ", "))

log_message("Samples: ", paste(unique(sample_values), collapse = ", "))

fib$fib_cluster_original <- factor(original_cluster, levels = expected_original_clusters)

merged_cluster <- original_cluster

merged_cluster[merged_cluster == "9"] <- "1"

expected_merged_clusters <- as.character(0:8)

if (!setequal(unique(merged_cluster), expected_merged_clusters)) {
    stop("Merged clusters are not exactly 0-8. Observed: ", paste(sort(unique(merged_cluster)), collapse = ", "))
}

fib$fib_cluster_merged <- factor(merged_cluster, levels = expected_merged_clusters)

expected_cluster1_n <- sum(original_cluster %in% c("1", "9"))

observed_cluster1_n <- sum(as.character(fib$fib_cluster_merged) == "1")

if (expected_cluster1_n != observed_cluster1_n) {
    stop("Cluster 9 -> 1 merge validation failed.")
}

log_message("Cluster 9 merged into cluster 1; merged cluster 1 cells: ", observed_cluster1_n)

subtype_levels <- c("PRG4+ lining fibroblasts", "APOD+ fibroblasts", "CD34+ sublining fibroblasts", "CXCL12+ sublining fibroblasts", 
    "DKK3+ sublining fibroblasts", "RSPO3+ fibroblasts", "POSTN+ fibroblasts")

merged_cluster_to_subtype <- c(`0` = "PRG4+ lining fibroblasts", `1` = "APOD+ fibroblasts", `2` = "CD34+ sublining fibroblasts", 
    `3` = "CXCL12+ sublining fibroblasts", `4` = "DKK3+ sublining fibroblasts", `5` = "PRG4+ lining fibroblasts", `6` = "CXCL12+ sublining fibroblasts", 
    `7` = "RSPO3+ fibroblasts", `8` = "POSTN+ fibroblasts")

final_subtype <- unname(merged_cluster_to_subtype[as.character(fib$fib_cluster_merged)])

if (anyNA(final_subtype)) {
    failed_clusters <- unique(as.character(fib$fib_cluster_merged)[is.na(final_subtype)])
    stop("No annotation was assigned to merged cluster(s): ", paste(failed_clusters, collapse = ", "))
}

fib$fibroblast_subtype_published <- factor(final_subtype, levels = subtype_levels)

if (any(table(fib$fibroblast_subtype_published) == 0)) {
    stop("At least one final subtype has zero cells.")
}

Idents(fib) <- "fibroblast_subtype_published"

annotation_mapping <- data.frame(original_cluster = expected_original_clusters, merged_cluster = c("0", "1", "2", "3", "4", 
    "5", "6", "7", "8", "1"), published_subtype = c("PRG4+ lining fibroblasts", "APOD+ fibroblasts", "CD34+ sublining fibroblasts", 
    "CXCL12+ sublining fibroblasts", "DKK3+ sublining fibroblasts", "PRG4+ lining fibroblasts", "CXCL12+ sublining fibroblasts", 
    "RSPO3+ fibroblasts", "POSTN+ fibroblasts", "APOD+ fibroblasts"), stringsAsFactors = FALSE)

annotation_mapping <- annotation_mapping[order(as.integer(annotation_mapping$original_cluster)), , drop = FALSE]

write.csv(annotation_mapping, file.path(table_dir, "01_cluster_to_published_subtype_mapping.csv"), row.names = FALSE, fileEncoding = "UTF-8")

subtype_counts <- as.data.frame(table(published_subtype = fib$fibroblast_subtype_published), stringsAsFactors = FALSE)

names(subtype_counts)[2] <- "n_cells"

subtype_counts$percent_of_fibroblasts <- subtype_counts$n_cells/sum(subtype_counts$n_cells) * 100

write.csv(subtype_counts, file.path(table_dir, "02_final_subtype_cell_counts_and_proportions.csv"), row.names = FALSE, fileEncoding = "UTF-8")

all_reductions <- Reductions(fib)

priority_reductions <- c("umap_fib_clean_corrected", "umap_fib_corrected", "umap_corrected", "umap_fib_clean", "umap_fib", 
    "umap")

umap_reduction <- priority_reductions[priority_reductions %in% all_reductions][1]

if (is.na(umap_reduction)) {
    umap_candidates <- all_reductions[grepl("umap", all_reductions, ignore.case = TRUE)]
    if (length(umap_candidates) == 0) {
        stop("No existing UMAP reduction was found. Available reductions: ", paste(all_reductions, collapse = ", "))
    }
    umap_reduction <- tail(umap_candidates, 1)
}

umap_embeddings <- Embeddings(fib, reduction = umap_reduction)

if (ncol(umap_embeddings) < 2) {
    stop("Selected UMAP reduction has fewer than two dimensions: ", umap_reduction)
}

if (!identical(rownames(umap_embeddings), colnames(fib))) {
    stop("UMAP embeddings and Seurat cell barcodes are not in identical order.")
}

log_message("Existing UMAP reduction used: ", umap_reduction)

log_message("No PCA, neighbors, clustering, or UMAP was recalculated.")

invisible(NULL)

invisible(NULL)

reference_palette <- c("#4E79A7", "#59A14F", "#F28E2B", "#E15759", "#B07AA1", "#76B7B2", "#EDC948", "#FF9DA7", "#9C755F")

cluster_palette <- setNames(reference_palette[seq_along(expected_merged_clusters)], expected_merged_clusters)

subtype_palette <- c(`PRG4+ lining fibroblasts` = "#4E79A7", `APOD+ fibroblasts` = "#59A14F", `CD34+ sublining fibroblasts` = "#F28E2B", 
    `CXCL12+ sublining fibroblasts` = "#E15759", `DKK3+ sublining fibroblasts` = "#B07AA1", `RSPO3+ fibroblasts` = "#FF9DA7", 
    `POSTN+ fibroblasts` = "#9C755F")

if (is.factor(fib[[sample_column]][, 1])) {
    sample_levels <- levels(droplevels(fib[[sample_column]][, 1]))
} else {
    sample_levels <- unique(sample_values)
}

fib$sample_id_for_umap <- factor(sample_values, levels = sample_levels)

if (length(sample_levels) <= length(reference_palette)) {
    sample_colors <- reference_palette[seq_along(sample_levels)]
} else {
    sample_colors <- grDevices::hcl.colors(length(sample_levels), palette = "Dark 3")
}

sample_palette <- setNames(sample_colors, sample_levels)

theme_reference_umap <- function(base_size = 12, legend_text_size = 10.5) {
    theme_classic(base_size = base_size, base_family = "Arial") + theme(plot.title = element_blank(), axis.title = element_text(size = 12, 
        colour = "black", face = "plain"), axis.text = element_blank(), axis.ticks = element_blank(), axis.line = element_line(colour = "black", 
        linewidth = 0.55), legend.position = "right", legend.title = element_text(size = 11, face = "plain", colour = "black"), 
        legend.text = element_text(size = legend_text_size, colour = "black"), legend.key = element_blank(), legend.background = element_blank(), 
        panel.background = element_rect(fill = "white", colour = NA), plot.background = element_rect(fill = "white", colour = NA), 
        plot.margin = margin(12, 12, 10, 10))
}

save_plot_pair <- function(plot_object, filename_without_extension, width = 7.2, height = 5.5, dpi = 600) {
    png_file <- paste0(filename_without_extension, ".png")
    pdf_file <- paste0(filename_without_extension, ".pdf")
    ggsave(filename = png_file, plot = plot_object, width = width, height = height, units = "in", dpi = dpi, bg = "white")
    if (capabilities("cairo")) {
        grDevices::cairo_pdf(filename = pdf_file, width = width, height = height, family = "Arial", onefile = TRUE)
    }
    else {
        grDevices::pdf(file = pdf_file, width = width, height = height, family = "sans", useDingbats = FALSE, onefile = TRUE)
    }
    print(plot_object)
    grDevices::dev.off()
    log_message("Saved: ", png_file)
    log_message("Saved: ", pdf_file)
}

p_cluster <- DimPlot(object = fib, reduction = umap_reduction, group.by = "fib_cluster_merged", cols = cluster_palette, pt.size = 0.25, 
    shuffle = TRUE, seed = 20260731, label = FALSE, raster = FALSE) + labs(x = "UMAP 1", y = "UMAP 2", colour = "Cluster") + 
    guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1))) + theme_reference_umap()

save_plot_pair(plot_object = p_cluster, filename_without_extension = file.path(figure_dir, "01_UMAP_merged_clusters_VECTOR"))

p_sample <- DimPlot(object = fib, reduction = umap_reduction, group.by = "sample_id_for_umap", cols = sample_palette, pt.size = 0.25, 
    shuffle = TRUE, seed = 20260731, label = FALSE, raster = FALSE) + labs(x = "UMAP 1", y = "UMAP 2", colour = "Sample") + 
    guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1))) + theme_reference_umap()

save_plot_pair(plot_object = p_sample, filename_without_extension = file.path(figure_dir, "02_UMAP_samples_VECTOR"))

p_annotation <- DimPlot(object = fib, reduction = umap_reduction, group.by = "fibroblast_subtype_published", cols = subtype_palette, 
    pt.size = 0.25, shuffle = TRUE, seed = 20260731, label = FALSE, raster = FALSE) + labs(x = "UMAP 1", y = "UMAP 2", colour = "Fibroblast subtype") + 
    guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1))) + theme_reference_umap(base_size = 12, legend_text_size = 9.5)

save_plot_pair(plot_object = p_annotation, filename_without_extension = file.path(figure_dir, "03_UMAP_published_fibroblast_subtypes_VECTOR"))

marker_panel <- data.frame(target_subtype = c(rep("PRG4+ lining fibroblasts", 4), rep("APOD+ fibroblasts", 4), rep("CD34+ sublining fibroblasts", 
    4), rep("CXCL12+ sublining fibroblasts", 4), rep("DKK3+ sublining fibroblasts", 4), rep("RSPO3+ fibroblasts", 3), rep("POSTN+ fibroblasts", 
    4)), gene = c("PRG4", "CLIC5", "HBEGF", "MMP3", "APOD", "CXCL14", "RARRES2", "IGFBP3", "CD34", "PI16", "DPP4", "COL15A1", 
    "CXCL12", "SFRP1", "CCL2", "ADAMTS1", "DKK3", "COMP", "PRELP", "OGN", "RSPO3", "APOC1", "APOE", "POSTN", "TNC", "COL6A1", 
    "LOXL2"), evidence_role = c(rep("published identity/reference marker", 20), "published identity marker", "current FindAllMarkers support; not used for naming", 
    "current FindAllMarkers support; not used for naming", "published identity marker", rep("published ECM-program support", 
        3)), stringsAsFactors = FALSE)

marker_panel$present_in_SCT <- marker_panel$gene %in% rownames(fib[["SCT"]])

missing_dotplot_genes <- marker_panel$gene[!marker_panel$present_in_SCT]

if (length(missing_dotplot_genes) > 0) {
    stop("Final DotPlot marker(s) absent from the SCT assay: ", paste(missing_dotplot_genes, collapse = ", "))
}

marker_genes <- marker_panel$gene

dot_calculation <- DotPlot(object = fib, features = marker_genes, assay = "SCT", group.by = "fibroblast_subtype_published", 
    scale = TRUE, col.min = -1, col.max = 2, dot.scale = 9, scale.by = "radius")

dot_data <- dot_calculation$data

required_dot_columns <- c("features.plot", "id", "pct.exp", "avg.exp", "avg.exp.scaled")

missing_dot_columns <- setdiff(required_dot_columns, colnames(dot_data))

if (length(missing_dot_columns) > 0) {
    stop("Unexpected Seurat DotPlot data structure. Missing column(s): ", paste(missing_dot_columns, collapse = ", "))
}

dot_data$gene <- factor(as.character(dot_data$features.plot), levels = rev(marker_genes))

dot_data$published_subtype <- factor(as.character(dot_data$id), levels = subtype_levels)

if (anyNA(dot_data$gene) || anyNA(dot_data$published_subtype)) {
    stop("DotPlot factor ordering failed.")
}

subtype_display_labels <- c(`PRG4+ lining fibroblasts` = "PRG4+ lining\nfibroblasts", `APOD+ fibroblasts` = "APOD+\nfibroblasts", 
    `CD34+ sublining fibroblasts` = "CD34+ sublining\nfibroblasts", `CXCL12+ sublining fibroblasts` = "CXCL12+ sublining\nfibroblasts", 
    `DKK3+ sublining fibroblasts` = "DKK3+ sublining\nfibroblasts", `RSPO3+ fibroblasts` = "RSPO3+\nfibroblasts", `POSTN+ fibroblasts` = "POSTN+\nfibroblasts")

p_dotplot <- ggplot(dot_data, aes(x = published_subtype, y = gene)) + geom_point(aes(size = pct.exp, colour = avg.exp.scaled)) + 
    scale_x_discrete(limits = subtype_levels, labels = subtype_display_labels, drop = FALSE) + scale_y_discrete(limits = rev(marker_genes), 
    drop = FALSE) + scale_size_continuous(name = "Per cent\nexpressed", limits = c(0, 100), breaks = c(25, 50, 75), range = c(0, 
    9)) + scale_colour_gradient2(name = "Average\nexpression", low = "#B2ABD2", mid = "#FFF7F3", high = "#E31A1C", midpoint = 0, 
    limits = c(-1, 2), breaks = c(-1, 0, 1, 2), oob = scales::squish) + labs(x = NULL, y = NULL) + guides(size = guide_legend(order = 1, 
    title.position = "top", override.aes = list(colour = "black")), colour = guide_colourbar(order = 2, title.position = "top", 
    title.hjust = 0.5, barheight = grid::unit(38, "mm"), barwidth = grid::unit(5.5, "mm"))) + theme_classic(base_size = 15, 
    base_family = "Arial") + theme(plot.title = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, 
    size = 13.5, colour = "black"), axis.text.y = element_text(size = 14.5, colour = "black"), axis.ticks = element_line(colour = "black", 
    linewidth = 0.45), axis.line = element_line(colour = "black", linewidth = 0.65), panel.grid = element_blank(), legend.position = "right", 
    legend.box = "vertical", legend.title = element_text(size = 13.5, face = "plain", colour = "black", lineheight = 1), 
    legend.text = element_text(size = 12.5, colour = "black"), legend.spacing.y = grid::unit(6, "mm"), legend.key = element_blank(), 
    legend.background = element_blank(), panel.background = element_rect(fill = "white", colour = NA), plot.background = element_rect(fill = "white", 
        colour = NA), plot.margin = margin(12, 12, 12, 12))

save_plot_pair(plot_object = p_dotplot, filename_without_extension = file.path(figure_dir, "04_DotPlot_final_published_fibroblast_subtypes"), 
    width = 9.2, height = 10.5, dpi = 600)

marker_match <- match(as.character(dot_data$gene), marker_panel$gene)

dotplot_source_data <- data.frame(published_subtype = as.character(dot_data$published_subtype), gene = as.character(dot_data$gene), 
    target_subtype = marker_panel$target_subtype[marker_match], evidence_role = marker_panel$evidence_role[marker_match], 
    percent_expressed = dot_data$pct.exp, average_expression = dot_data$avg.exp, scaled_average_expression = dot_data$avg.exp.scaled, 
    stringsAsFactors = FALSE)

write.csv(dotplot_source_data, file.path(table_dir, "03_final_subtype_DotPlot_source_data.csv"), row.names = FALSE, fileEncoding = "UTF-8")

log_message("Final annotation DotPlot and source-data CSV saved.")

fib[["sample_id_for_umap"]] <- NULL

output_rds <- file.path(rds_dir, "GSE216651_Fibroblasts_final_7_published_subtypes.rds")

saveRDS(fib, file = output_rds, compress = TRUE)

log_message("Annotated RDS saved: ", output_rds)

log_message("Final subtype counts:")

for (i in seq_len(nrow(subtype_counts))) {
    log_message("  ", subtype_counts$published_subtype[i], ": ", subtype_counts$n_cells[i], " cells (", sprintf("%.2f", subtype_counts$percent_of_fibroblasts[i]), 
        "%)")
}

log_message("Completed: ", format(Sys.time()))

cat("\nAnalysis completed successfully.\n")

cat("Output directory: ", output_dir, "\n", sep = "")

cat("Annotated RDS: ", output_rds, "\n", sep = "")

cat("UMAP reduction reused: ", umap_reduction, "\n", sep = "")

cat("No PCA, neighbors, clustering, or UMAP was recalculated.\n")

