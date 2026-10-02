args <- commandArgs(trailingOnly = FALSE)

script_file <- sub("^--file=", "", grep("^--file=", args, value = TRUE)[1])

code_root <- dirname(dirname(normalizePath(script_file, mustWork = TRUE)))

source(file.path(code_root, "00_configuration.R"))

suppressPackageStartupMessages({
    library(Seurat)
    library(Matrix)
})

set.seed(20260712)

project_dir <- normalizePath(file.path(input_root, "GSE216651_singlecell"))

raw_dir <- file.path(project_dir, "GSE216651_RAW", "extracted")

metadata_dir <- file.path(project_dir, "00_metadata")

output_dir <- file.path(project_dir, "01_scRNA_preprocess", "01_raw_seurat_objects")

log_dir <- file.path(project_dir, "summary", "logs")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

manifest <- read.delim(file.path(metadata_dir, "GSE216651_sample_manifest.tsv"), stringsAsFactors = FALSE, check.names = FALSE)

manifest <- manifest[manifest$modality == "scRNA", , drop = FALSE]

if (nrow(manifest) != 6L) stop("Expected exactly six scRNA-seq samples; found ", nrow(manifest))

if (!identical(sort(unique(manifest$group)), c("Control", "OA"))) stop("Expected Control and OA groups")

target_sample <- Sys.getenv("TARGET_SAMPLE", unset = "")

if (nzchar(target_sample)) {
    manifest <- manifest[manifest$sample_id == target_sample, , drop = FALSE]
    if (nrow(manifest) != 1L) 
        stop("Unknown TARGET_SAMPLE: ", target_sample)
}

summary_rows <- vector("list", nrow(manifest))

for (i in seq_len(nrow(manifest))) {
    x <- manifest[i, ]
    prefix <- file.path(raw_dir, paste0(x$GSM, "_", x$sample_id))
    paths <- c(mtx = paste0(prefix, "_matrix.mtx.gz"), features = paste0(prefix, "_features.tsv.gz"), barcodes = paste0(prefix, 
        "_barcodes.tsv.gz"))
    if (!all(file.exists(paths))) 
        stop("Missing 10x file(s) for ", x$sample_id)
    counts <- ReadMtx(mtx = paths[["mtx"]], features = paths[["features"]], cells = paths[["barcodes"]], feature.column = 2, 
        cell.column = 1, unique.features = TRUE)
    obj <- CreateSeuratObject(counts = counts, project = x$sample_id, min.cells = 0, min.features = 0, names.field = 1, names.delim = "_")
    obj <- RenameCells(obj, add.cell.id = x$sample_id)
    obj$sample_id <- x$sample_id
    obj$donor_id <- x$donor_id
    obj$Group <- x$group
    obj$GSM <- x$GSM
    obj$modality <- x$modality
    obj$source_tissue <- x$source_tissue
    obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
    if (ncol(obj) == 0L) 
        stop("No barcodes loaded for ", x$sample_id)
    if (!all(startsWith(colnames(obj), paste0(x$sample_id, "_")))) 
        stop("Barcode prefix failed for ", x$sample_id)
    saveRDS(obj, file.path(output_dir, paste0(x$sample_id, "_raw_unfiltered_seurat.rds")), compress = FALSE)
    cell_metadata <- data.frame(cell_id = colnames(obj), obj[[]][, c("orig.ident", "nCount_RNA", "nFeature_RNA", "percent.mt", 
        "sample_id", "donor_id", "Group", "GSM", "modality", "source_tissue"), drop = FALSE], check.names = FALSE)
    write.table(cell_metadata, gzfile(file.path(output_dir, paste0(x$sample_id, "_raw_unfiltered_cell_metadata.tsv.gz"))), 
        sep = "\t", quote = FALSE, row.names = FALSE)
    summary_rows[[i]] <- data.frame(sample_id = x$sample_id, GSM = x$GSM, donor_id = x$donor_id, Group = x$group, modality = x$modality, 
        raw_cells = ncol(obj), raw_features = nrow(obj), median_nCount_RNA = median(obj$nCount_RNA), median_nFeature_RNA = median(obj$nFeature_RNA), 
        median_percent_mt = median(obj$percent.mt), stringsAsFactors = FALSE)
    rm(counts, obj, cell_metadata)
    gc()
}

summary_table <- do.call(rbind, summary_rows)

write.table(summary_table, file.path(output_dir, if (nzchar(target_sample)) paste0(target_sample, "_raw_seurat_object_summary.tsv") else "scRNA_raw_seurat_object_summary.tsv"), 
    sep = "\t", quote = FALSE, row.names = FALSE)

if (!nzchar(target_sample)) {
    write.table(manifest, file.path(output_dir, "scRNA_analysis_sample_manifest.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
}

writeLines(capture.output(sessionInfo()), file.path(log_dir, "01_create_scRNA_seurat_objects_sessionInfo.txt"))

message("Created unfiltered Seurat objects for ", nrow(manifest), " scRNA-seq samples; total cells = ", sum(summary_table$raw_cells))

