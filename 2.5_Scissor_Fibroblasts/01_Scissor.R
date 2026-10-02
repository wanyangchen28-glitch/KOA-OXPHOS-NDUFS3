args <- commandArgs(trailingOnly = FALSE)

script_file <- sub("^--file=", "", grep("^--file=", args, value = TRUE)[1])

code_root <- dirname(dirname(normalizePath(script_file, mustWork = TRUE)))

source(file.path(code_root, "00_configuration.R"))

suppressPackageStartupMessages({
    library(Seurat)
    library(Matrix)
    library(dplyr)
    library(Scissor)
    library(preprocessCore)
})

SAMPLE_FIELD <- "sample_id"

CELLTYPE_FIELD <- "major_cell_type"

RANDOM_SEED <- 20260730L

proportional_stratified_sample <- function(meta, target_n, seed) {
    stopifnot(target_n > 0L, target_n <= nrow(meta))
    allocation <- meta %>% count(.data[[SAMPLE_FIELD]], .data[[CELLTYPE_FIELD]], name = "stratum_n") %>% mutate(exact_n = stratum_n/sum(stratum_n) * 
        target_n, take_n = floor(exact_n), remainder = exact_n - take_n)
    cells_left <- target_n - sum(allocation$take_n)
    if (cells_left > 0L) {
        eligible <- which(allocation$take_n < allocation$stratum_n)
        add_order <- eligible[order(allocation$remainder[eligible], decreasing = TRUE)]
        allocation$take_n[add_order[seq_len(cells_left)]] <- allocation$take_n[add_order[seq_len(cells_left)]] + 1L
    }
    set.seed(seed)
    selected <- vector("list", nrow(allocation))
    for (i in seq_len(nrow(allocation))) {
        in_stratum <- rownames(meta)[meta[[SAMPLE_FIELD]] == allocation[[SAMPLE_FIELD]][i] & meta[[CELLTYPE_FIELD]] == allocation[[CELLTYPE_FIELD]][i]]
        selected[[i]] <- if (allocation$take_n[i] > 0L) {
            sample(in_stratum, size = allocation$take_n[i], replace = FALSE)
        }
        else {
            character(0)
        }
    }
    selected_cells <- unlist(selected, use.names = FALSE)
    if (length(selected_cells) != target_n) {
        stop("Proportional sampling returned ", length(selected_cells), " cells instead of ", target_n, ".")
    }
    list(cells = selected_cells, allocation = allocation)
}

run_scissor_v5 <- function(bulk_dataset, sc_expression, network_sparse, phenotype, alpha_search, cutoff, seed) {
    common <- intersect(rownames(bulk_dataset), rownames(sc_expression))
    if (length(common) == 0L) {
        stop("No common genes between bulk and single-cell matrices.")
    }
    bulk_common <- as.matrix(bulk_dataset[common, , drop = FALSE])
    sc_common_sparse <- sc_expression[common, , drop = FALSE]
    detected_sc <- Matrix::rowSums(sc_common_sparse != 0) > 0
    finite_bulk <- apply(bulk_common, 1, function(z) all(is.finite(z)))
    keep <- detected_sc & finite_bulk
    common <- common[keep]
    bulk_common <- bulk_common[keep, , drop = FALSE]
    sc_common_sparse <- sc_common_sparse[keep, , drop = FALSE]
    message("Common detected genes used by Scissor: ", length(common))
    message("Converting sampled single-cell expression to a dense matrix.")
    sc_dense <- as.matrix(sc_common_sparse)
    dataset0 <- cbind(bulk_common, sc_dense)
    rm(sc_dense, sc_common_sparse, bulk_common)
    invisible(gc())
    message("Performing joint quantile normalization.")
    dataset1 <- preprocessCore::normalize.quantiles(dataset0, copy = FALSE)
    rownames(dataset1) <- common
    colnames(dataset1) <- colnames(dataset0)
    rm(dataset0)
    invisible(gc())
    bulk_n <- ncol(bulk_dataset)
    expression_bulk <- dataset1[, seq_len(bulk_n), drop = FALSE]
    expression_cell <- dataset1[, (bulk_n + 1L):ncol(dataset1), drop = FALSE]
    message("Calculating bulk-cell Pearson correlation matrix.")
    X <- stats::cor(expression_bulk, expression_cell, method = "pearson")
    quality_check <- stats::quantile(X, na.rm = TRUE)
    rm(dataset1, expression_bulk, expression_cell)
    invisible(gc())
    if (anyNA(X) || any(!is.finite(X))) {
        stop("The bulk-cell correlation matrix contains non-finite values.")
    }
    if (quality_check[[3]] < 0.01) {
        warning("Median bulk-cell correlation is below 0.01.")
    }
    message("Converting the sampled SNN graph to the Scissor network matrix.")
    network_sparse@x[] <- 1
    diag(network_sparse) <- 0
    network_sparse <- Matrix::drop0(network_sparse)
    network <- as.matrix(network_sparse)
    rm(network_sparse)
    invisible(gc())
    if (!identical(colnames(X), rownames(network))) {
        network <- network[colnames(X), colnames(X), drop = FALSE]
    }
    Y <- as.numeric(phenotype)
    if (!identical(sort(unique(Y)), c(0, 1))) {
        stop("Phenotype must contain exactly 0 (Control) and 1 (OA).")
    }
    chosen <- NULL
    set.seed(seed)
    for (alpha_value in alpha_search) {
        message("Fitting Scissor model: alpha = ", alpha_value)
        fit_cv <- Scissor::APML1(X, Y, family = "binomial", penalty = "Net", alpha = alpha_value, Omega = network, nlambda = 100, 
            nfolds = min(10, nrow(X)))
        fit_final <- Scissor::APML1(X, Y, family = "binomial", penalty = "Net", alpha = alpha_value, Omega = network, lambda = fit_cv$lambda.min)
        coefficients <- as.numeric(fit_final$Beta[2:(ncol(X) + 1L)])
        names(coefficients) <- colnames(X)
        positive_cells <- names(coefficients)[coefficients > 0]
        negative_cells <- names(coefficients)[coefficients < 0]
        selected_fraction <- (length(positive_cells) + length(negative_cells))/length(coefficients)
        message("  Scissor+ = ", length(positive_cells), "; Scissor- = ", length(negative_cells), "; selected = ", sprintf("%.3f%%", 
            selected_fraction * 100))
        chosen <- list(para = list(alpha = alpha_value, lambda = fit_cv$lambda.min, family = "binomial", cutoff = cutoff), 
            Coefs = coefficients, Scissor_pos = positive_cells, Scissor_neg = negative_cells, selected_fraction = selected_fraction, 
            quality_check = quality_check, common_genes = common)
        if (selected_fraction < cutoff) 
            break
    }
    rm(X, network)
    invisible(gc())
    chosen
}

bulk_dir <- file.path(input_root, "bulk_discovery", "00_preprocess")

bulk <- as.matrix(read.csv(file.path(bulk_dir, "06_ComBat_HGNC_maxIQR_expression.csv"), row.names = 1, check.names = FALSE))

meta_bulk <- read.delim(file.path(bulk_dir, "04_sample_info.tsv"))

meta_bulk <- meta_bulk[match(colnames(bulk), meta_bulk$GSM), , drop = FALSE]

phenotype <- as.integer(meta_bulk$Group == "OA")

seu <- readRDS(file.path(input_root, "singlecell", "GSE216651_scRNA_PC30_res0.6_without_clusters16_17.rds"))

sampling <- proportional_stratified_sample(seu[[]], 8000L, RANDOM_SEED)

seu <- subset(seu, cells = sampling$cells)

DefaultAssay(seu) <- "SCT"

seu <- FindNeighbors(seu, reduction = "pca", dims = 1:30, k.param = 20, graph.name = c("Scissor_nn", "Scissor_snn"), verbose = FALSE)

result <- run_scissor_v5(bulk, GetAssayData(seu, assay = "SCT", layer = "data"), seu@graphs[["Scissor_snn"]], phenotype, 
    c(0.005, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9), 0.2, RANDOM_SEED)

meta <- seu[[]]

meta$cell_id <- rownames(meta)

meta$Scissor_coefficient <- result$Coefs[meta$cell_id]

meta$Scissor_raw_class <- "Background"

meta$Scissor_raw_class[meta$cell_id %in% result$Scissor_pos] <- "Scissor+"

meta$Scissor_raw_class[meta$cell_id %in% result$Scissor_neg] <- "Scissor-"

outdir <- file.path(output_root, "Scissor")

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

write.csv(meta, file.path(outdir, "cell_labels.csv"), row.names = FALSE)

write.csv(sampling$allocation, file.path(outdir, "sampling_allocation.csv"), row.names = FALSE)

saveRDS(result, file.path(outdir, "Scissor_result.rds"))

