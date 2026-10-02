args <- commandArgs(trailingOnly = FALSE)

script_file <- sub("^--file=", "", grep("^--file=", args, value = TRUE)[1])

code_root <- dirname(dirname(normalizePath(script_file, mustWork = TRUE)))

source(file.path(code_root, "00_configuration.R"))

suppressPackageStartupMessages({
    library(Seurat)
    library(SingleCellExperiment)
    library(scDblFinder)
    library(BiocParallel)
})

project_dir <- normalizePath(file.path(input_root, "GSE216651_singlecell"))

input_dir <- file.path(project_dir, "01_scRNA_preprocess", "02_original_QC_filtered")

output_dir <- file.path(project_dir, "01_scRNA_preprocess", "03_scDblFinder_filtered")

metadata_dir <- file.path(project_dir, "00_metadata")

log_dir <- file.path(project_dir, "summary", "logs")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

manifest <- read.delim(file.path(metadata_dir, "GSE216651_sample_manifest.tsv"), stringsAsFactors = FALSE, check.names = FALSE)

manifest <- manifest[manifest$modality == "scRNA", , drop = FALSE]

if (nrow(manifest) != 6L) stop("Expected exactly six scRNA-seq samples; found ", nrow(manifest))

target_sample <- Sys.getenv("TARGET_SAMPLE", unset = "")

if (nzchar(target_sample)) {
    manifest <- manifest[manifest$sample_id == target_sample, , drop = FALSE]
    if (nrow(manifest) != 1L) 
        stop("Unknown TARGET_SAMPLE: ", target_sample)
}

params <- list(scDblFinder_version = as.character(packageVersion("scDblFinder")), clusters = NULL, samples = NULL, dbr = NULL, 
    dbr_per_1k = 0.008, nfeatures = 1352, dims = 20, score = "xgb", threshold = TRUE, processing = "default", BPPARAM = "SerialParam")

summary_rows <- vector("list", nrow(manifest))

for (i in seq_len(nrow(manifest))) {
    x <- manifest[i, ]
    set.seed(20260712L + i)
    input_rds <- file.path(input_dir, paste0(x$sample_id, "_original_QC_filtered_seurat.rds"))
    if (!file.exists(input_rds)) 
        stop("Missing QC-filtered Seurat object: ", input_rds)
    obj <- readRDS(input_rds)
    sce <- as.SingleCellExperiment(obj, assay = "RNA")
    if (!identical(colnames(obj), colnames(sce))) 
        stop("Cell order changed during SCE conversion: ", x$sample_id)
    sce <- scDblFinder(sce, clusters = NULL, samples = NULL, dbr = NULL, dbr.per1k = 0.008, nfeatures = 1352, dims = 20, 
        score = "xgb", threshold = TRUE, processing = "default", BPPARAM = SerialParam(progressbar = FALSE), verbose = TRUE)
    if (!all(c("scDblFinder.score", "scDblFinder.class") %in% colnames(colData(sce)))) 
        stop("scDblFinder output missing for ", x$sample_id)
    if (!identical(colnames(obj), colnames(sce))) 
        stop("Cell order changed after scDblFinder: ", x$sample_id)
    obj$scDblFinder.score <- colData(sce)$scDblFinder.score
    obj$scDblFinder.class <- as.character(colData(sce)$scDblFinder.class)
    if (anyNA(obj$scDblFinder.class)) 
        stop("Missing scDblFinder class: ", x$sample_id)
    if (!all(obj$scDblFinder.class %in% c("singlet", "doublet"))) 
        stop("Unexpected scDblFinder class: ", x$sample_id)
    decisions <- data.frame(cell_id = colnames(obj), obj[[]][, c("orig.ident", "nCount_RNA", "nFeature_RNA", "percent.mt", 
        "sample_id", "donor_id", "Group", "GSM", "scDblFinder.score", "scDblFinder.class"), drop = FALSE], keep_after_scDblFinder = obj$scDblFinder.class == 
        "singlet", check.names = FALSE)
    write.table(decisions, gzfile(file.path(output_dir, paste0(x$sample_id, "_scDblFinder_cell_decisions.tsv.gz"))), sep = "\t", 
        quote = FALSE, row.names = FALSE)
    singlet_obj <- subset(obj, subset = scDblFinder.class == "singlet")
    if (ncol(singlet_obj) < 3L) 
        stop("Fewer than three singlets remain: ", x$sample_id)
    saveRDS(singlet_obj, file.path(output_dir, paste0(x$sample_id, "_QC_scDblFinder_singlets_seurat.rds")), compress = FALSE)
    class_counts <- table(obj$scDblFinder.class)
    summary_rows[[i]] <- data.frame(sample_id = x$sample_id, GSM = x$GSM, donor_id = x$donor_id, Group = x$group, cells_before_scDblFinder = ncol(obj), 
        singlets = unname(class_counts["singlet"]), doublets = unname(class_counts["doublet"]), doublet_percent = 100 * unname(class_counts["doublet"])/ncol(obj), 
        median_singlet_score = median(obj$scDblFinder.score[obj$scDblFinder.class == "singlet"]), median_doublet_score = median(obj$scDblFinder.score[obj$scDblFinder.class == 
            "doublet"]), random_seed = 20260712L + i, stringsAsFactors = FALSE)
    rm(obj, sce, singlet_obj, decisions)
    gc()
}

summary_table <- do.call(rbind, summary_rows)

summary_name <- if (nzchar(target_sample)) paste0(target_sample, "_scDblFinder_summary.tsv") else "scRNA_scDblFinder_summary.tsv"

write.table(summary_table, file.path(output_dir, summary_name), sep = "\t", quote = FALSE, row.names = FALSE)

parameter_table <- data.frame(parameter = names(params), value = vapply(params, function(v) if (is.null(v)) "NULL" else paste(v, 
    collapse = ","), character(1)), stringsAsFactors = FALSE)

write.table(parameter_table, file.path(output_dir, "scDblFinder_parameter_manifest.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

writeLines(capture.output(sessionInfo()), file.path(log_dir, "04_scDblFinder_per_donor_sessionInfo.txt"))

message("Completed scDblFinder for ", nrow(manifest), " scRNA-seq sample(s).")

