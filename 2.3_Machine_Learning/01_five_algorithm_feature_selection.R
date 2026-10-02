args <- commandArgs(trailingOnly = FALSE)

script_file <- sub("^--file=", "", grep("^--file=", args, value = TRUE)[1])

code_root <- dirname(dirname(normalizePath(script_file, mustWork = TRUE)))

source(file.path(code_root, "00_configuration.R"))

options(stringsAsFactors = FALSE, width = 160)

project_dir <- file.path(input_root, "bulk_training")

outdir <- file.path(output_root, "feature_selection")

source_dir <- file.path(outdir, "source_data")

code_dir <- file.path(outdir, "code")

log_dir <- file.path(project_dir, "summary", "logs")

local_lib <- file.path(project_dir, "R_library")

candidate_file <- file.path(input_root, "bulk_discovery/10_OXPHOS_bulk_GSEA/08_OXPHOS_DEG_WGCNA_intersection_candidates.csv")

train_expression_file <- file.path(project_dir, "00_preprocess", "06_expanded_training_ComBat_HGNC_expression.csv")

train_metadata_file <- file.path(project_dir, "00_preprocess", "group_metadata.csv")

seed_lasso <- 11L

seed_enet <- 12L

seed_rf <- 3L

seed_boruta <- 1L

seed_svmrfe <- 123L

seed_roc_cv <- 2026L

rf_ntree <- 2000L

boruta_max_runs <- 500L

cv_folds <- 5L

lasso_lambda_rule <- "lambda.1se"

expected_group_counts <- c(Control = 27L, OA = 39L)

condition_colours <- c(Control = "#2D6AA0", OA = "#C94B4B")

required_packages <- c("dplyr", "tidyr", "tibble", "ggplot2", "glmnet", "randomForest", "Boruta", "caret", "e1071", "pROC", 
    "ragg")

missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing_packages) > 0L) stop("Missing R packages: ", paste(missing_packages, collapse = ", "))

suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(tibble)
    library(ggplot2)
    library(glmnet)
    library(randomForest)
    library(Boruta)
    library(caret)
    library(e1071)
    library(pROC)
    library(ragg)
})

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)

dir.create(code_dir, recursive = TRUE, showWarnings = FALSE)

dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

dir.create(local_lib, recursive = TRUE, showWarnings = FALSE)

.libPaths(c(local_lib, .libPaths()))

theme_article <- function(base_size = 8.2, base_family = "Arial") {
    theme_classic(base_size = base_size, base_family = base_family) + theme(axis.line = element_line(linewidth = 0.45, colour = "black"), 
        axis.ticks = element_line(linewidth = 0.4, colour = "black"), axis.ticks.length = grid::unit(1.8, "pt"), axis.title = element_text(size = base_size, 
            colour = "black"), axis.text = element_text(size = base_size - 0.7, colour = "black"), legend.title = element_blank(), 
        legend.text = element_text(size = base_size - 0.5, colour = "black"), plot.title = element_text(size = base_size + 
            1.1, face = "bold", hjust = 0), plot.subtitle = element_text(size = base_size - 0.5, colour = "#404040", hjust = 0), 
        panel.grid = element_blank(), plot.margin = ggplot2::margin(5.5, 7, 5, 5.5))
}

save_gg <- function(plot, prefix, width_mm, height_mm, dpi = 600) {
    width_in <- width_mm/25.4
    height_in <- height_mm/25.4
    grDevices::cairo_pdf(paste0(prefix, ".pdf"), width = width_in, height = height_in, family = "Arial")
    print(plot)
    dev.off()
    ragg::agg_png(paste0(prefix, ".png"), width = width_in, height = height_in, units = "in", res = dpi)
    print(plot)
    dev.off()
}

save_base <- function(draw_fun, prefix, width_mm, height_mm, dpi = 600) {
    width_in <- width_mm/25.4
    height_in <- height_mm/25.4
    grDevices::cairo_pdf(paste0(prefix, ".pdf"), width = width_in, height = height_in, family = "Arial")
    draw_fun()
    dev.off()
    ragg::agg_png(paste0(prefix, ".png"), width = width_in, height = height_in, units = "in", res = dpi)
    draw_fun()
    dev.off()
}

read_expression <- function(path) {
    x <- read.csv(path, check.names = FALSE, row.names = 1, stringsAsFactors = FALSE)
    x <- as.matrix(x)
    storage.mode(x) <- "numeric"
    x
}

if (!all(file.exists(c(candidate_file, train_expression_file, train_metadata_file)))) {
    missing_input <- c(candidate_file, train_expression_file, train_metadata_file)
    stop("Missing input: ", paste(missing_input[!file.exists(missing_input)], collapse = "; "))
}

candidates <- read.csv(candidate_file, check.names = FALSE, stringsAsFactors = FALSE)

if (!"gene" %in% names(candidates) && "Gene" %in% names(candidates)) names(candidates)[names(candidates) == "Gene"] <- "gene"

if (!"gene" %in% names(candidates)) stop("OXPHOS candidate file must contain a gene column.")

candidates$Gene <- as.character(candidates$gene)

candidates$strict_candidate <- TRUE

candidates <- candidates %>% distinct(Gene, .keep_all = TRUE)

candidate_genes <- c("NDUFB7", "NDUFA3", "NDUFA8", "ATP6V1A", "NDUFS3", "COX5B")

candidates <- candidates %>% filter(Gene %in% candidate_genes)

if (!identical(sort(candidates$Gene), sort(candidate_genes))) stop("OXPHOS turquoise candidate definition does not contain exactly the expected six genes.")

candidates <- candidates[match(candidate_genes, candidates$Gene), , drop = FALSE]

train_expr <- read_expression(train_expression_file)

train_meta <- read.csv(train_metadata_file, check.names = FALSE, stringsAsFactors = FALSE) %>% filter(Group %in% c("Control", 
    "OA")) %>% select(Sample, Group)

train_meta$Group <- factor(train_meta$Group, levels = c("Control", "OA"))

train_meta <- train_meta[match(colnames(train_expr), train_meta$Sample), , drop = FALSE]

if (anyNA(train_meta$Group) || !identical(train_meta$Sample, colnames(train_expr))) stop("Training expression/metadata alignment failed.")

if (!identical(as.integer(table(train_meta$Group)[names(expected_group_counts)]), as.integer(expected_group_counts))) stop("Unexpected expanded-training group counts.")

if (!all(candidate_genes %in% rownames(train_expr))) stop("Training matrix missing candidate genes: ", paste(setdiff(candidate_genes, 
    rownames(train_expr)), collapse = ", "))

x_df <- as.data.frame(t(train_expr[candidate_genes, train_meta$Sample, drop = FALSE]), check.names = FALSE)

x <- as.matrix(x_df)

storage.mode(x) <- "numeric"

y_factor <- factor(train_meta$Group, levels = c("Control", "OA"))

y <- as.integer(y_factor == "OA")

if (nrow(x) < 3L || ncol(x) < 1L || min(table(y_factor)) < cv_folds) stop("Insufficient samples for configured 5-fold CV.")

training_input <- data.frame(Sample = train_meta$Sample, Group = y_factor, x_df, check.names = FALSE)

write.csv(training_input, file.path(source_dir, "01_expanded_training_ML_input_sample_by_gene.csv"), row.names = FALSE)

write.csv(candidates, file.path(source_dir, "02_frozen_OXPHOS_six_candidate_definition.csv"), row.names = FALSE)

set.seed(seed_lasso)

cv_lasso <- glmnet::cv.glmnet(x, y, alpha = 1, family = "binomial", nfolds = cv_folds, type.measure = "deviance")

coef_lasso <- as.matrix(coef(cv_lasso, s = lasso_lambda_rule))

lasso_df <- data.frame(Gene = rownames(coef_lasso), Coefficient = as.numeric(coef_lasso[, 1]), stringsAsFactors = FALSE) %>% 
    filter(Gene != "(Intercept)", Coefficient != 0) %>% arrange(Coefficient)

lasso_genes <- lasso_df$Gene

write.csv(data.frame(Log_lambda = log(cv_lasso$lambda), Mean_deviance = cv_lasso$cvm, SE_deviance = cv_lasso$cvsd, N_nonzero = cv_lasso$nzero), 
    file.path(source_dir, "03_LASSO_CV_results.csv"), row.names = FALSE)

write.csv(lasso_df, file.path(source_dir, "04_LASSO_selected_coefficients.csv"), row.names = FALSE)

save_base(function() {
    par(family = "Arial", cex.axis = 0.9, cex.lab = 0.95, cex.main = 0.9)
    plot(cv_lasso, main = "LASSO cross-validation")
}, file.path(outdir, "01_LASSO_CV_curve"), 150, 115)

set.seed(seed_enet)

cv_enet <- glmnet::cv.glmnet(x, y, alpha = 0.5, family = "binomial", nfolds = cv_folds, type.measure = "deviance")

coef_enet <- as.matrix(coef(cv_enet, s = "lambda.min"))

enet_df <- data.frame(Gene = rownames(coef_enet), Coefficient = as.numeric(coef_enet[, 1]), stringsAsFactors = FALSE) %>% 
    filter(Gene != "(Intercept)", Coefficient != 0) %>% arrange(Coefficient)

enet_genes <- enet_df$Gene

write.csv(data.frame(Log_lambda = log(cv_enet$lambda), Mean_deviance = cv_enet$cvm, SE_deviance = cv_enet$cvsd, N_nonzero = cv_enet$nzero), 
    file.path(source_dir, "05_Elastic_Net_CV_results.csv"), row.names = FALSE)

write.csv(enet_df, file.path(source_dir, "06_Elastic_Net_selected_coefficients.csv"), row.names = FALSE)

p_enet <- ggplot(enet_df, aes(x = reorder(Gene, Coefficient), y = Coefficient)) + geom_col(fill = "#5A5A5A", width = 0.72) + 
    geom_hline(yintercept = 0, linewidth = 0.45) + coord_flip() + labs(title = "Elastic Net selected gene coefficients", 
    x = NULL, y = "Coefficient") + theme_article() + theme(axis.text.y = element_text(face = "bold"))

save_gg(p_enet, file.path(outdir, "04_Elastic_Net_selected_gene_coefficients"), 150, 105)

set.seed(seed_rf)

rf_fit <- randomForest::randomForest(x = x_df, y = y_factor, importance = TRUE, ntree = rf_ntree, proximity = TRUE)

rf_importance <- as.data.frame(randomForest::importance(rf_fit, scale = FALSE)) %>% rownames_to_column("Gene") %>% arrange(desc(MeanDecreaseGini))

rf_genes <- rf_importance$Gene[seq_len(min(15L, nrow(rf_importance)))]

write.csv(rf_importance, file.path(source_dir, "07_Random_Forest_feature_importance.csv"), row.names = FALSE)

p_rf <- ggplot(rf_importance, aes(x = reorder(Gene, MeanDecreaseGini), y = MeanDecreaseGini)) + geom_col(fill = "#4C9F70", 
    width = 0.72) + coord_flip() + labs(title = "Random Forest variable importance", x = NULL, y = "MeanDecreaseGini") + 
    theme_article() + theme(axis.text.y = element_text(face = "bold"))

save_gg(p_rf, file.path(outdir, "05_Random_Forest_variable_importance"), 150, 105)

boruta_input <- cbind(data.frame(Group = y_factor), x_df)

set.seed(seed_boruta)

boruta_fit <- Boruta::Boruta(Group ~ ., data = boruta_input, doTrace = 0, maxRuns = boruta_max_runs)

boruta_stats <- Boruta::attStats(boruta_fit) %>% as.data.frame() %>% rownames_to_column("Gene") %>% arrange(desc(meanImp))

boruta_genes <- boruta_stats %>% filter(decision == "Confirmed") %>% pull(Gene)

write.csv(boruta_stats, file.path(source_dir, "09_Boruta_feature_statistics.csv"), row.names = FALSE)

save_base(function() {
    par(family = "Arial", mar = c(10.5, 4.5, 3.5, 1.5), cex.axis = 0.7, cex.lab = 0.95, cex.main = 0.9)
    plot(boruta_fit, las = 2, cex.axis = 0.7, xlab = "", main = "Boruta feature importance")
}, file.path(outdir, "07_Boruta_feature_importance"), 170, 120)

set.seed(seed_svmrfe)

svm_ctrl <- caret::rfeControl(functions = caret::caretFuncs, method = "cv", number = cv_folds, returnResamp = "final")

svm_profile <- caret::rfe(x = x_df, y = y_factor, sizes = seq_len(ncol(x_df)), rfeControl = svm_ctrl, method = "svmLinear")

svm_result <- svm_profile$results %>% arrange(Variables)

svm_best_n <- svm_profile$optsize

svm_genes <- caret::predictors(svm_profile)

write.csv(svm_result, file.path(source_dir, "10_SVM_RFE_cross_validation_results.csv"), row.names = FALSE)

write.csv(data.frame(Gene = svm_genes), file.path(source_dir, "11_SVM_RFE_selected_features.csv"), row.names = FALSE)

best_svm <- svm_result %>% filter(Variables == svm_best_n) %>% slice(1)

p_svm <- ggplot(svm_result, aes(x = Variables, y = Accuracy)) + geom_line(colour = "#8FC9E2", linewidth = 0.85) + geom_point(colour = "#8FC9E2", 
    size = 1.8) + geom_point(data = best_svm, colour = "#D6604D", size = 2.8) + geom_text(data = best_svm, aes(label = paste0("Best: ", 
    Variables, " feature", ifelse(Variables > 1, "s", ""), "\nAccuracy = ", sprintf("%.3f", Accuracy))), vjust = 1.5, size = 2.45, 
    colour = "#202020") + scale_x_continuous(breaks = svm_result$Variables) + labs(title = "SVM-RFE cross-validation", x = "Number of selected features", 
    y = "Accuracy") + theme_article()

save_gg(p_svm, file.path(outdir, "09_SVM_RFE_cross_validation"), 150, 110)

ml_list <- list(LASSO = unique(lasso_genes), ElasticNet = unique(enet_genes), RF = unique(rf_genes), Boruta = unique(boruta_genes), 
    SVMRFE = unique(svm_genes))

membership <- expand.grid(Gene = candidate_genes, Method = names(ml_list), stringsAsFactors = FALSE) %>% mutate(Selected = mapply(function(g, 
    m) g %in% ml_list[[m]], Gene, Method))

selection_summary <- membership %>% group_by(Gene) %>% summarise(Methods_selected = sum(Selected), Selected_by = paste(Method[Selected], 
    collapse = "; "), .groups = "drop") %>% arrange(desc(Methods_selected), Gene)

strict_intersection <- Reduce(intersect, ml_list)

write.csv(membership, file.path(source_dir, "12_ML_selection_membership.csv"), row.names = FALSE)

write.csv(selection_summary, file.path(source_dir, "13_ML_gene_selection_frequency.csv"), row.names = FALSE)

write.csv(data.frame(Gene = strict_intersection), file.path(source_dir, "14_ML_strict_all_five_method_intersection.csv"), 
    row.names = FALSE)
