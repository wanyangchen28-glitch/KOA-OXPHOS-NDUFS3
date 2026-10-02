args <- commandArgs(trailingOnly = FALSE)

script_file <- sub("^--file=", "", grep("^--file=", args, value = TRUE)[1])

code_root <- dirname(dirname(normalizePath(script_file, mustWork = TRUE)))

source(file.path(code_root, "00_configuration.R"))

suppressPackageStartupMessages(library(Seurat))

project_dir <- normalizePath(file.path(input_root, "GSE216651_singlecell"))

input_dir <- file.path(project_dir, "01_scRNA_preprocess", "01_raw_seurat_objects")

output_dir <- file.path(project_dir, "01_scRNA_preprocess", "02_original_QC_filtered")

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

summary_rows <- vector("list", nrow(manifest))

for (i in seq_len(nrow(manifest))) {
    x <- manifest[i, ]
    input_rds <- file.path(input_dir, paste0(x$sample_id, "_raw_unfiltered_seurat.rds"))
    if (!file.exists(input_rds)) 
        stop("Missing first-step Seurat object: ", input_rds)
    obj <- readRDS(input_rds)
    meta <- obj[[]]
    fail_low_feature <- meta$nFeature_RNA < 200
    fail_high_feature <- meta$nFeature_RNA > 6000
    fail_high_mt <- meta$percent.mt > 25
    keep <- !(fail_low_feature | fail_high_feature | fail_high_mt)
    cell_qc <- data.frame(cell_id = rownames(meta), meta[, c("orig.ident", "nCount_RNA", "nFeature_RNA", "percent.mt", "sample_id", 
        "donor_id", "Group", "GSM", "modality", "source_tissue"), drop = FALSE], fail_nFeature_lt_200 = fail_low_feature, 
        fail_nFeature_gt_6000 = fail_high_feature, fail_percent_mt_gt_25 = fail_high_mt, keep_original_QC = keep, check.names = FALSE)
    write.table(cell_qc, gzfile(file.path(output_dir, paste0(x$sample_id, "_original_QC_cell_decisions.tsv.gz"))), sep = "\t", 
        quote = FALSE, row.names = FALSE)
    filtered_obj <- subset(obj, cells = rownames(meta)[keep])
    if (ncol(filtered_obj) < 3L) 
        stop("Fewer than three cells remain after QC for ", x$sample_id)
    saveRDS(filtered_obj, file.path(output_dir, paste0(x$sample_id, "_original_QC_filtered_seurat.rds")), compress = FALSE)
    summary_rows[[i]] <- data.frame(sample_id = x$sample_id, GSM = x$GSM, donor_id = x$donor_id, Group = x$group, raw_cells = ncol(obj), 
        excluded_nFeature_lt_200 = sum(fail_low_feature), excluded_nFeature_gt_6000 = sum(fail_high_feature), excluded_percent_mt_gt_25 = sum(fail_high_mt), 
        excluded_any_original_QC = sum(!keep), retained_original_QC = sum(keep), retained_percent = 100 * mean(keep), median_nCount_RNA_after_QC = median(filtered_obj$nCount_RNA), 
        median_nFeature_RNA_after_QC = median(filtered_obj$nFeature_RNA), median_percent_mt_after_QC = median(filtered_obj$percent.mt), 
        stringsAsFactors = FALSE)
    rm(obj, filtered_obj, meta, cell_qc)
    gc()
}

summary_table <- do.call(rbind, summary_rows)

summary_name <- if (nzchar(target_sample)) paste0(target_sample, "_original_QC_filter_summary.tsv") else "scRNA_original_QC_filter_summary.tsv"

write.table(summary_table, file.path(output_dir, summary_name), sep = "\t", quote = FALSE, row.names = FALSE)

if (!nzchar(target_sample)) {
    invisible(NULL)
}

writeLines(capture.output(sessionInfo()), file.path(log_dir, "02_filter_scRNA_original_QC_sessionInfo.txt"))

message("Applied original QC thresholds to ", nrow(manifest), " scRNA-seq sample(s).")

