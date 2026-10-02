args <- commandArgs(trailingOnly = FALSE)

script_file <- sub("^--file=", "", grep("^--file=", args, value = TRUE)[1])

code_root <- dirname(dirname(normalizePath(script_file, mustWork = TRUE)))

source(file.path(code_root, "00_configuration.R"))

suppressPackageStartupMessages({
    library(Seurat)
    library(ggplot2)
})

project_dir <- normalizePath(file.path(input_root, "GSE216651_singlecell"))

input_dir <- file.path(project_dir, "01_scRNA_preprocess", "04_cellcycle_SCT_per_donor")

output_dir <- file.path(project_dir, "01_scRNA_preprocess", "05_SCT_integration_PCA_unselected")

metadata_dir <- file.path(project_dir, "00_metadata")

log_dir <- file.path(project_dir, "summary", "logs")

figure_dir <- file.path(output_dir, "figures")

source_dir <- file.path(figure_dir, "source_data")

dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)

dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

manifest <- read.delim(file.path(metadata_dir, "GSE216651_sample_manifest.tsv"), stringsAsFactors = FALSE, check.names = FALSE)

manifest <- manifest[manifest$modality == "scRNA", , drop = FALSE]

objects <- lapply(manifest$sample_id, function(sample_id) {
    f <- file.path(input_dir, paste0(sample_id, "_cellcycle_SCT_seurat.rds"))
    if (!file.exists(f)) 
        stop("Missing SCT object: ", f)
    readRDS(f)
})

names(objects) <- manifest$sample_id

features <- SelectIntegrationFeatures(object.list = objects, nfeatures = 2000)

objects <- PrepSCTIntegration(object.list = objects, anchor.features = features, verbose = TRUE)

anchors <- FindIntegrationAnchors(object.list = objects, normalization.method = "SCT", anchor.features = features, reduction = "cca", 
    dims = 1:30, verbose = TRUE)

integrated <- IntegrateData(anchorset = anchors, normalization.method = "SCT", dims = 1:30, new.assay.name = "integrated", 
    verbose = TRUE)

DefaultAssay(integrated) <- "integrated"

set.seed(20260712L)

integrated <- RunPCA(integrated, assay = "integrated", npcs = 50, seed.use = 20260712L, verbose = TRUE)

pc_stdev <- integrated[["pca"]]@stdev

pc_table <- data.frame(PC = seq_along(pc_stdev), stdev = pc_stdev, variance = pc_stdev^2, percent_variance = 100 * pc_stdev^2/sum(pc_stdev^2), 
    cumulative_percent_variance = 100 * cumsum(pc_stdev^2)/sum(pc_stdev^2), stringsAsFactors = FALSE)

write.table(features, file.path(output_dir, "SCT_integration_features_2000.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)

write.table(pc_table, file.path(output_dir, "PCA_variance_50_components.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

write.table(data.frame(parameter = c("SelectIntegrationFeatures_nfeatures", "FindIntegrationAnchors_normalization.method", 
    "FindIntegrationAnchors_reduction", "FindIntegrationAnchors_dims", "IntegrateData_normalization.method", "IntegrateData_dims", 
    "RunPCA_npcs", "RunPCA_seed.use"), value = c("2000", "SCT", "cca", "1:30", "SCT", "1:30", "50", "20260712"), stringsAsFactors = FALSE), 
    file.path(output_dir, "SCT_integration_PCA_parameter_manifest.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

saveRDS(integrated, file.path(output_dir, "GSE216651_scRNA_SCT_integrated_PCA50_unselected.rds"), compress = FALSE)

p <- ggplot(pc_table, aes(x = PC, y = stdev)) + geom_line(linewidth = 0.5, colour = "#1F77B4") + geom_point(size = 1.5, colour = "#1F77B4") + 
    scale_x_continuous(breaks = seq(0, 50, 5), limits = c(1, 50)) + labs(title = "PCA elbow plot after SCT anchor integration", 
    subtitle = "50 PCs calculated; no downstream PC cutoff has been selected.", x = "Principal component", y = "Standard deviation") + 
    theme_classic(base_size = 9, base_family = "Arial") + theme(axis.line = element_line(linewidth = 0.4), axis.ticks = element_line(linewidth = 0.35), 
    plot.title = element_text(face = "bold", size = 11), plot.subtitle = element_text(size = 8, colour = "#4D4D4D"))

write.table(pc_table, file.path(source_dir, "PCA_elbow_source_data.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

width_mm <- 100

height_mm <- 75

width_in <- width_mm/25.4

height_in <- height_mm/25.4

stub <- file.path(figure_dir, "Fig_scRNA_SCT_integrated_PCA50_elbow_unselected")

grDevices::cairo_pdf(paste0(stub, ".pdf"), width = width_in, height = height_in, family = "Arial")

print(p)

grDevices::dev.off()

ragg::agg_tiff(paste0(stub, ".tiff"), width = width_in, height = height_in, units = "in", res = 600, compression = "lzw")

print(p)

grDevices::dev.off()

writeLines(capture.output(sessionInfo()), file.path(log_dir, "07_SCT_anchor_integration_PCA_unselected_sessionInfo.txt"))

message("Completed SCT integration and PCA(50); no PC cutoff selected.")

