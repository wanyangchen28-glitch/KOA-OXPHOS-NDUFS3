args <- commandArgs(trailingOnly = FALSE)

script_file <- sub("^--file=", "", grep("^--file=", args, value = TRUE)[1])

code_root <- dirname(dirname(normalizePath(script_file, mustWork = TRUE)))

source(file.path(code_root, "00_configuration.R"))

options(stringsAsFactors = FALSE, width = 160)

project_root <- input_root

dataset_dir <- file.path(project_root, "GSE89408_validation")

count_file <- file.path(dataset_dir, "GSE89408_RAW", "GSE89408_GEO_count_matrix_rename.txt.gz")

preprocess_dir <- file.path(dataset_dir, "00_preprocess")

dir.create(preprocess_dir, recursive = TRUE, showWarnings = FALSE)

alias_map <- c(CECR1 = "ADA2")

expected_group_counts <- c(Control = 28L, OA = 22L)

required_packages <- c("data.table", "edgeR", "pROC")

missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing_packages) > 0) stop("Missing required packages: ", paste(missing_packages, collapse = ", "))

if (!file.exists(count_file)) stop("Count matrix not found: ", count_file)

message("Reading GSE89408 supplementary matrix: ", basename(count_file))

count_dt <- data.table::fread(count_file, data.table = FALSE, check.names = FALSE)

if (ncol(count_dt) < 2) stop("The supplementary matrix has no sample columns.")

input_gene <- trimws(as.character(count_dt[[1]]))

all_samples <- colnames(count_dt)[-1]

selected_samples <- all_samples[grepl("^(normal_tissue|OA_tissue)_[0-9]+$", all_samples)]

if (length(selected_samples) != sum(expected_group_counts)) {
    stop("Expected 50 OA/Normal samples but selected ", length(selected_samples), ".")
}

metadata <- data.frame(Sample = selected_samples, Group = ifelse(grepl("^OA_tissue_", selected_samples), "OA", "Control"), 
    Dataset = "GSE89408", MatrixSample = selected_samples, stringsAsFactors = FALSE)

metadata$Group <- factor(metadata$Group, levels = c("Control", "OA"))

if (!identical(as.integer(table(metadata$Group)[c("Control", "OA")]), as.integer(expected_group_counts))) {
    stop("Unexpected group counts after sample selection.")
}

raw_matrix <- as.matrix(count_dt[, selected_samples, drop = FALSE])

storage.mode(raw_matrix) <- "numeric"

rownames(raw_matrix) <- input_gene

if (anyNA(raw_matrix) || any(!is.finite(raw_matrix))) stop("Input count matrix contains missing or non-finite values.")

if (any(raw_matrix < 0)) stop("Input count matrix contains negative values.")

keep_nonzero <- rowSums(raw_matrix > 0) > 0

raw_matrix <- raw_matrix[keep_nonzero, , drop = FALSE]

input_gene <- input_gene[keep_nonzero]

mapped_gene <- input_gene

mapped_gene[input_gene %in% names(alias_map)] <- unname(alias_map[input_gene[input_gene %in% names(alias_map)]])

if (any(!nzchar(mapped_gene))) stop("Empty gene symbols remain after mapping.")

dge <- edgeR::DGEList(counts = raw_matrix)

dge <- edgeR::calcNormFactors(dge, method = "TMM")

log_cpm <- edgeR::cpm(dge, normalized.lib.sizes = TRUE, log = TRUE, prior.count = 1)

log_cpm_sum <- rowsum(log_cpm, group = mapped_gene, reorder = FALSE)

mapped_counts <- table(mapped_gene)

log_cpm_hgnc <- sweep(log_cpm_sum, 1, as.numeric(mapped_counts[rownames(log_cpm_sum)]), "/")

rownames(log_cpm_hgnc) <- rownames(log_cpm_sum)

write.csv(metadata, file.path(preprocess_dir, "01_GSE89408_OA_Normal_metadata.csv"), row.names = FALSE)

write.csv(log_cpm_hgnc, file.path(preprocess_dir, "03_GSE89408_TMM_log2CPM_HGNC_expression.csv"), row.names = TRUE)

