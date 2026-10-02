args <- commandArgs(trailingOnly = FALSE)

script_file <- sub("^--file=", "", grep("^--file=", args, value = TRUE)[1])

code_root <- dirname(dirname(normalizePath(script_file, mustWork = TRUE)))

source(file.path(code_root, "00_configuration.R"))

suppressPackageStartupMessages(library(pROC))

dataset_dir <- file.path(input_root, "GSE89408_validation", "00_preprocess")

outdir <- file.path(output_root, "GSE89408_validation")

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

expr <- as.matrix(read.csv(file.path(dataset_dir, "03_GSE89408_TMM_log2CPM_HGNC_expression.csv"), row.names = 1, check.names = FALSE))

meta <- read.csv(file.path(dataset_dir, "01_GSE89408_OA_Normal_metadata.csv"))

meta <- meta[match(colnames(expr), meta$MatrixSample), , drop = FALSE]

stopifnot(!anyNA(meta$Group), sum(meta$Group == "OA") == 22, sum(meta$Group == "Control") == 28)

cf <- read.csv(file.path(output_root, "two_gene_model", "09_ATP6V1A_NDUFS3_frozen_coefficients_and_training_scaling.csv"))

coefficients <- setNames(cf$Frozen_coefficient, cf$Gene_or_Term)

genes <- c("ATP6V1A", "NDUFS3")

z <- scale(t(expr[genes, , drop = FALSE]))

score <- as.numeric(z[, genes, drop = FALSE] %*% coefficients[genes] + coefficients["(Intercept)"])

y <- as.integer(meta$Group == "OA")

ro <- roc(y, score, levels = c(0, 1), direction = "<", quiet = TRUE)

ci <- as.numeric(ci.auc(ro, method = "delong"))

metrics <- data.frame(Cohort = "GSE89408", N = length(y), Control = sum(y == 0), KOA = sum(y == 1), AUC = as.numeric(auc(ro)), 
    CI_low = ci[1], CI_high = ci[3])

write.csv(metrics, file.path(outdir, "ROC_metrics.csv"), row.names = FALSE)

write.csv(data.frame(Sample = meta$Sample, Group = meta$Group, Score = score), file.path(outdir, "predictions.csv"), row.names = FALSE)

write.csv(coords(ro, "all", ret = c("threshold", "specificity", "sensitivity"), transpose = FALSE), file.path(outdir, "ROC_coordinates.csv"), 
    row.names = FALSE)

