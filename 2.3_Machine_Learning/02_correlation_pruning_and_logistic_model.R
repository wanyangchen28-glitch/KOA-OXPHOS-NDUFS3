args <- commandArgs(trailingOnly = FALSE)

script_file <- sub("^--file=", "", grep("^--file=", args, value = TRUE)[1])

code_root <- dirname(dirname(normalizePath(script_file, mustWork = TRUE)))

source(file.path(code_root, "00_configuration.R"))

suppressPackageStartupMessages({
    library(dplyr)
    library(tibble)
    library(pROC)
})

project_dir <- file.path(input_root, "bulk_training")

source_dir <- file.path(output_root, "two_gene_model")

dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)

consensus_file <- file.path(output_root, "feature_selection", "source_data", "14_ML_strict_all_five_method_intersection.csv")

three_gene_consensus <- read.csv(consensus_file)$Gene

redundancy_threshold <- 0.8

train_expr <- as.matrix(read.csv(file.path(project_dir, "00_preprocess", "06_expanded_training_ComBat_HGNC_expression.csv"), 
    row.names = 1, check.names = FALSE))

train_meta <- read.csv(file.path(project_dir, "00_preprocess", "group_metadata.csv"))

train_meta <- train_meta[match(colnames(train_expr), train_meta$Sample), , drop = FALSE]

train_meta$Group <- factor(train_meta$Group, levels = c("Control", "OA"))

stopifnot(identical(train_meta$Sample, colnames(train_expr)), !anyNA(train_meta$Group))

candidate_file <- file.path(output_root, "OXPHOS_GSEA", "candidate_genes.csv")

candidate_stats <- read.csv(candidate_file)

candidate_stats$Gene <- candidate_stats$gene

expected_retained_genes <- c("ATP6V1A", "NDUFS3")

zscore_fit_apply <- function(train_x, test_x = NULL) {
    train_x <- as.matrix(train_x)
    if (is.null(test_x)) 
        test_x <- train_x
    test_x <- as.matrix(test_x)
    center <- colMeans(train_x, na.rm = TRUE)
    scale_value <- apply(train_x, 2, sd, na.rm = TRUE)
    scale_value[!is.finite(scale_value) | scale_value == 0] <- 1
    list(train = sweep(sweep(train_x, 2, center, "-"), 2, scale_value, "/"), test = sweep(sweep(test_x, 2, center, "-"), 
        2, scale_value, "/"), center = center, scale = scale_value)
}

roc_summary <- function(labels, probabilities, cohort_label) {
    roc_obj <- pROC::roc(labels, probabilities, levels = c(0, 1), direction = "<", quiet = TRUE)
    ci <- as.numeric(pROC::ci.auc(roc_obj, method = "delong"))
    coords <- pROC::coords(roc_obj, x = "all", ret = c("threshold", "specificity", "sensitivity"), transpose = FALSE)
    curve <- data.frame(Model = cohort_label, Threshold = coords$threshold, Specificity = coords$specificity, Sensitivity = coords$sensitivity, 
        FPR = 1 - coords$specificity, TPR = coords$sensitivity, stringsAsFactors = FALSE) %>% arrange(FPR, TPR)
    metrics <- data.frame(Cohort = cohort_label, AUC = as.numeric(pROC::auc(roc_obj)), AUC_CI_low = ci[1], AUC_CI_high = ci[3], 
        N = length(labels), Control = sum(labels == 0), OA = sum(labels == 1), stringsAsFactors = FALSE)
    list(curve = curve, metrics = metrics)
}

train_three <- t(train_expr[three_gene_consensus, train_meta$Sample, drop = FALSE])

cor_matrix <- cor(train_three, method = "pearson", use = "pairwise.complete.obs")

write.csv(cor_matrix, file.path(source_dir, "02_training_only_three_gene_Pearson_correlation_matrix.csv"), row.names = TRUE)

distance_matrix <- as.dist(1 - abs(cor_matrix))

hc <- hclust(distance_matrix, method = "average")

cluster_id <- cutree(hc, h = 1 - redundancy_threshold)

cluster_df <- data.frame(Gene = names(cluster_id), Cluster = as.integer(cluster_id), stringsAsFactors = FALSE)

cluster_df <- cluster_df %>% left_join(candidate_stats, by = "Gene")

representatives <- cluster_df %>% group_by(Cluster) %>% arrange(bulk_FDR, .by_group = TRUE) %>% slice(1) %>% ungroup() %>% 
    pull(Gene)

representatives <- representatives[match(three_gene_consensus[three_gene_consensus %in% representatives], representatives)]

representatives <- unique(representatives)

if (!setequal(representatives, expected_retained_genes)) {
    stop("Training-only correlation pruning did not yield the expected ATP6V1A-NDUFS3 pair. Retained: ", paste(representatives, 
        collapse = ", "))
}

cluster_df$Retained_representative <- cluster_df$Gene %in% representatives

write.csv(cluster_df, file.path(source_dir, "03_training_only_redundancy_cluster_assignment.csv"), row.names = FALSE)

train_x <- t(train_expr[expected_retained_genes, train_meta$Sample, drop = FALSE])

train_x <- as.data.frame(train_x, check.names = FALSE)

train_y <- as.integer(train_meta$Group == "OA")

scaled_train <- zscore_fit_apply(train_x)

scaled_train_df <- as.data.frame(scaled_train$train, check.names = FALSE)

fit_data <- data.frame(Label = train_y, scaled_train_df, check.names = FALSE)

fit_formula <- as.formula(paste("Label ~", paste(expected_retained_genes, collapse = " + ")))

signature_fit <- stats::glm(fit_formula, data = fit_data, family = stats::binomial())

fit_coef <- as.data.frame(summary(signature_fit)$coefficients) %>% rownames_to_column("Term")

names(fit_coef) <- c("Term", "Estimate", "Std_Error", "z_value", "P_value")

fit_coef$Model <- "ATP6V1A_NDUFS3_compact_logistic"

fit_coef <- fit_coef %>% select(Model, Term, Estimate, Std_Error, z_value, P_value)

write.csv(fit_coef, file.path(source_dir, "05_ATP6V1A_NDUFS3_frozen_logistic_coefficients.csv"), row.names = FALSE)

train_probability <- as.numeric(stats::predict(signature_fit, newdata = scaled_train_df, type = "response"))

train_apparent <- roc_summary(train_y, train_probability, "66-sample training (apparent)")

frozen_scaling <- data.frame(Gene_or_Term = c("(Intercept)", expected_retained_genes), Training_center = c(NA, scaled_train$center[expected_retained_genes]), 
    Training_scale = c(NA, scaled_train$scale[expected_retained_genes]), Frozen_coefficient = c(coef(signature_fit)["(Intercept)"], 
        coef(signature_fit)[expected_retained_genes]), stringsAsFactors = FALSE)

write.csv(frozen_scaling, file.path(source_dir, "09_ATP6V1A_NDUFS3_frozen_coefficients_and_training_scaling.csv"), row.names = FALSE)

write.csv(train_apparent$metrics, file.path(source_dir, "training_ROC_metrics.csv"), row.names = FALSE)

write.csv(train_apparent$curve, file.path(source_dir, "training_ROC_coordinates.csv"), row.names = FALSE)

write.csv(data.frame(Sample = train_meta$Sample, Group = train_meta$Group, Probability = train_probability), file.path(source_dir, 
    "training_predictions.csv"), row.names = FALSE)

saveRDS(signature_fit, file.path(source_dir, "logistic_model.rds"))

