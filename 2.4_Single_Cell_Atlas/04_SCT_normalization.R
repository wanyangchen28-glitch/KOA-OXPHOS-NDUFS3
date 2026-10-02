args <- commandArgs(trailingOnly = FALSE)

script_file <- sub("^--file=", "", grep("^--file=", args, value = TRUE)[1])

code_root <- dirname(dirname(normalizePath(script_file, mustWork = TRUE)))

source(file.path(code_root, "00_configuration.R"))

suppressPackageStartupMessages(library(Seurat))

options(future.globals.maxSize = 2 * 1024^3)

project_dir <- normalizePath(file.path(input_root, "GSE216651_singlecell"))

input_dir <- file.path(project_dir, "01_scRNA_preprocess", "03_scDblFinder_filtered")

output_dir <- file.path(project_dir, "01_scRNA_preprocess", "04_cellcycle_SCT_per_donor")

metadata_dir <- file.path(project_dir, "00_metadata")

log_dir <- file.path(project_dir, "summary", "logs")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

manifest <- read.delim(file.path(metadata_dir, "GSE216651_sample_manifest.tsv"), stringsAsFactors = FALSE, check.names = FALSE)

manifest <- manifest[manifest$modality == "scRNA", , drop = FALSE]

target_sample <- Sys.getenv("TARGET_SAMPLE", unset = "")

if (nzchar(target_sample)) {
    manifest <- manifest[manifest$sample_id == target_sample, , drop = FALSE]
    if (nrow(manifest) != 1L) 
        stop("Unknown TARGET_SAMPLE: ", target_sample)
}

cc_s <- cc.genes.updated.2019$s.genes

cc_g2m <- cc.genes.updated.2019$g2m.genes

summary_rows <- vector("list", nrow(manifest))

for (i in seq_len(nrow(manifest))) {
    x <- manifest[i, ]
    input_rds <- file.path(input_dir, paste0(x$sample_id, "_QC_scDblFinder_singlets_seurat.rds"))
    if (!file.exists(input_rds)) 
        stop("Missing singlet Seurat object: ", input_rds)
    obj <- readRDS(input_rds)
    DefaultAssay(obj) <- "RNA"
    obj <- NormalizeData(obj, assay = "RNA", normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
    obj <- CellCycleScoring(obj, s.features = cc_s, g2m.features = cc_g2m, set.ident = FALSE)
    set.seed(20260712L)
    obj <- SCTransform(obj, assay = "RNA", new.assay.name = "SCT", ncells = 5000, variable.features.n = 3000, vars.to.regress = c("S.Score", 
        "G2M.Score"), vst.flavor = "v2", seed.use = 20260712L, return.only.var.genes = TRUE, verbose = TRUE)
    if (!all(c("S.Score", "G2M.Score", "Phase") %in% colnames(obj[[]]))) 
        stop("Cell-cycle metadata missing: ", x$sample_id)
    if (!"SCT" %in% names(obj@assays)) 
        stop("SCT assay missing: ", x$sample_id)
    phase_tab <- table(factor(obj$Phase, levels = c("G1", "S", "G2M")))
    summary_rows[[i]] <- data.frame(sample_id = x$sample_id, GSM = x$GSM, donor_id = x$donor_id, Group = x$group, singlet_cells = ncol(obj), 
        G1_cells = unname(phase_tab["G1"]), S_cells = unname(phase_tab["S"]), G2M_cells = unname(phase_tab["G2M"]), median_S_score = median(obj$S.Score), 
        median_G2M_score = median(obj$G2M.Score), SCT_variable_features = length(VariableFeatures(obj, assay = "SCT")), stringsAsFactors = FALSE)
    saveRDS(obj, file.path(output_dir, paste0(x$sample_id, "_cellcycle_SCT_seurat.rds")), compress = FALSE)
    rm(obj)
    gc()
}

summary_table <- do.call(rbind, summary_rows)

summary_name <- if (nzchar(target_sample)) paste0(target_sample, "_cellcycle_SCT_summary.tsv") else "scRNA_cellcycle_SCT_summary.tsv"

write.table(summary_table, file.path(output_dir, summary_name), sep = "\t", quote = FALSE, row.names = FALSE)

write.table(data.frame(parameter = c("temporary_NormalizeData", "cell_cycle_gene_set", "CellCycleScoring_set.ident", "SCTransform_assay", 
    "SCTransform_vars.to.regress", "SCTransform_ncells", "SCTransform_variable.features.n", "SCTransform_vst.flavor", "SCTransform_seed.use", 
    "SCTransform_return.only.var.genes"), value = c("LogNormalize; scale.factor=10000", "cc.genes.updated.2019", "FALSE", 
    "RNA", "S.Score,G2M.Score", "5000", "3000", "v2", "20260712", "TRUE"), stringsAsFactors = FALSE), file.path(output_dir, 
    "cellcycle_SCT_parameter_manifest.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

writeLines(capture.output(sessionInfo()), file.path(log_dir, "06_cellcycle_SCT_per_donor_sessionInfo.txt"))

message("Completed cell-cycle scoring and SCTransform for ", nrow(manifest), " sample(s).")

