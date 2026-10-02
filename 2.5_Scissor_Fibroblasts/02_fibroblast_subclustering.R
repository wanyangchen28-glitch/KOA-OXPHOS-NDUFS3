args <- commandArgs(trailingOnly = FALSE)

script_file <- sub("^--file=", "", grep("^--file=", args, value = TRUE)[1])

code_root <- dirname(dirname(normalizePath(script_file, mustWork = TRUE)))

source(file.path(code_root, "00_configuration.R"))

suppressPackageStartupMessages(library(Seurat))

obj <- readRDS(file.path(input_root, "singlecell", "GSE216651_scRNA_PC30_res0.6_without_clusters16_17.rds"))

fibro <- subset(obj, cells = colnames(obj)[obj$major_cell_type == "Fibroblasts"])

DefaultAssay(fibro) <- "integrated"

features <- VariableFeatures(fibro[["integrated"]])

if (!length(features)) features <- rownames(fibro[["integrated"]])

set.seed(20260730)

fibro <- RunPCA(fibro, assay = "integrated", features = features, npcs = 50, reduction.name = "fibro_pca50", reduction.key = "fPC_", 
    seed.use = 20260730)

fibro <- FindNeighbors(fibro, reduction = "fibro_pca50", dims = 1:20, k.param = 20, graph.name = c("fibro_PC20_nn", "fibro_PC20_snn"))

fibro <- RunUMAP(fibro, reduction = "fibro_pca50", dims = 1:20, n.neighbors = 30, min.dist = 0.3, metric = "cosine", umap.method = "uwot", 
    reduction.name = "fibro_umap_PC20", reduction.key = "fUMAP20_", seed.use = 20260730)

fibro <- FindClusters(fibro, graph.name = "fibro_PC20_snn", resolution = 0.5, algorithm = 1, random.seed = 20260730, cluster.name = "fibro_res_0_5")

outdir <- file.path(output_root, "fibroblast_clustering")

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

saveRDS(fibro, file.path(outdir, "fibroblast_PC20_res0.5.rds"))

