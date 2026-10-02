args <- commandArgs(trailingOnly = FALSE)

script_file <- sub("^--file=", "", grep("^--file=", args, value = TRUE)[1])

code_root <- dirname(dirname(normalizePath(script_file, mustWork = TRUE)))

source(file.path(code_root, "00_configuration.R"))

INPUT_RDS <- file.path(input_root, "singlecell/GSE216651_scRNA_PC30_res0.6_without_clusters16_17.rds")

OUT_DIR <- file.path(output_root, "36_NDUFS3_OA_CellChat_exact_template")

RDS_DIR <- file.path(OUT_DIR, "01_rds")

TABLE_DIR <- file.path(OUT_DIR, "02_tables")

TEMPLATE_FIG_DIR <- file.path(OUT_DIR, "03_template_figures")

OPTIONAL_FIG_DIR <- file.path(OUT_DIR, "04_optional_single_figures")

LOG_DIR <- file.path(OUT_DIR, "05_logs")

for (x in c(RDS_DIR, TABLE_DIR, TEMPLATE_FIG_DIR, OPTIONAL_FIG_DIR, LOG_DIR)) {
    dir.create(x, recursive = TRUE, showWarnings = FALSE)
}

ASSAY_NAME <- "SCT"

LAYER_NAME <- "data"

GROUP_FIELD <- "Group"

MAJOR_FIELD <- "major_cell_type"

OA_LABEL <- "OA"

TARGET_GENE <- "NDUFS3"

MIN_CELLS <- 10

WORKERS <- 10

RUN_PATTERN_ANALYSIS <- TRUE

N_PATTERNS <- 2

PATHWAY_SHOW <- "CCL"

CELLCHAT_LEVELS <- c("NDUFS3+ Fibroblasts", "NDUFS3- Fibroblasts", "Endothelial", "Mural", "Myeloid", "Lymphoid", "Mast")

FOCUS_GROUPS <- c("NDUFS3+ Fibroblasts", "NDUFS3- Fibroblasts")

FOCUS_COLORS <- c(`NDUFS3+ Fibroblasts` = "#D84A43", `NDUFS3- Fibroblasts` = "#4C78A8")

options(stringsAsFactors = FALSE)

set.seed(123)

required_packages <- c("Seurat", "SeuratObject", "CellChat", "future", "ggplot2", "RColorBrewer", "ComplexHeatmap", "circlize")

missing_packages <- required_packages[!vapply(required_packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]

if (length(missing_packages) > 0) {
    stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
    library(Seurat)
    library(SeuratObject)
    library(CellChat)
    library(future)
    library(ggplot2)
    library(RColorBrewer)
    library(ComplexHeatmap)
    library(circlize)
})

write_csv <- function(x, path) {
    utils::write.csv(x, path, row.names = FALSE, fileEncoding = "UTF-8")
}

safe_name <- function(x) {
    gsub("[^A-Za-z0-9]+", "_", x)
}

save_base_both <- function(stem, width, height, draw_fun, folder = TEMPLATE_FIG_DIR) {
    grDevices::cairo_pdf(file.path(folder, paste0(stem, ".pdf")), width, height)
    tryCatch(draw_fun(), finally = grDevices::dev.off())
    grDevices::png(file.path(folder, paste0(stem, ".png")), width = width, height = height, units = "in", res = 320, bg = "white")
    tryCatch(draw_fun(), finally = grDevices::dev.off())
}

save_gg_both <- function(plot_object, stem, width, height, folder = TEMPLATE_FIG_DIR) {
    ggplot2::ggsave(file.path(folder, paste0(stem, ".pdf")), plot_object, width = width, height = height, units = "in", device = grDevices::cairo_pdf, 
        bg = "white")
    ggplot2::ggsave(file.path(folder, paste0(stem, ".png")), plot_object, width = width, height = height, units = "in", dpi = 320, 
        bg = "white")
}

save_heatmap_both <- function(heatmap_object, stem, width, height) {
    save_base_both(stem, width, height, function() {
        ComplexHeatmap::draw(heatmap_object, heatmap_legend_side = "right", annotation_legend_side = "right")
    })
}

if (!file.exists(INPUT_RDS)) stop("Input RDS not found: ", INPUT_RDS)

object <- readRDS(INPUT_RDS)

if (!all(c(GROUP_FIELD, MAJOR_FIELD) %in% colnames(object@meta.data))) {
    stop("Group or major_cell_type metadata is missing.")
}

data.input.all <- SeuratObject::LayerData(object[[ASSAY_NAME]], layer = LAYER_NAME)

gene.index <- which(toupper(rownames(data.input.all)) == TARGET_GENE)

if (length(gene.index) != 1) stop("NDUFS3 was not uniquely found in SCT/data.")

oa.cells <- rownames(object@meta.data)[as.character(object@meta.data[[GROUP_FIELD]]) == OA_LABEL]

oa.meta <- object@meta.data[oa.cells, , drop = FALSE]

oa.meta$CellChat_group <- as.character(oa.meta[[MAJOR_FIELD]])

fib.index <- which(as.character(oa.meta[[MAJOR_FIELD]]) == "Fibroblasts")

fib.cells <- rownames(oa.meta)[fib.index]

ndufs3.value <- as.numeric(data.input.all[gene.index, fib.cells])

oa.meta$CellChat_group[fib.index[ndufs3.value > 0]] <- "NDUFS3+ Fibroblasts"

oa.meta$CellChat_group[fib.index[ndufs3.value == 0]] <- "NDUFS3- Fibroblasts"

oa.meta$CellChat_group <- factor(oa.meta$CellChat_group, levels = CELLCHAT_LEVELS)

if (anyNA(oa.meta$CellChat_group)) stop("Unexpected CellChat identity was found.")

oa.data.input <- data.input.all[, oa.cells, drop = FALSE]

if (!identical(colnames(oa.data.input), rownames(oa.meta))) {
    stop("Expression columns and metadata rows are not aligned.")
}

cell.counts <- as.data.frame(table(oa.meta$CellChat_group))

colnames(cell.counts) <- c("CellChat_group", "n_cells")

write_csv(cell.counts, file.path(TABLE_DIR, "01_CellChat_group_counts.csv"))

rm(object, data.input.all)

invisible(gc())

oa.cellchat <- CellChat::createCellChat(object = oa.data.input)

oa.cellchat <- CellChat::addMeta(oa.cellchat, meta = oa.meta)

oa.cellchat <- CellChat::setIdent(oa.cellchat, ident.use = "CellChat_group")

groupSize <- as.numeric(table(oa.cellchat@idents))

names(groupSize) <- levels(oa.cellchat@idents)

oa.cellchat@DB <- CellChat::CellChatDB.human

oa.cellchat <- CellChat::subsetData(oa.cellchat, features = NULL)

future::plan("multisession", workers = WORKERS)

oa.cellchat <- CellChat::identifyOverExpressedGenes(oa.cellchat)

oa.cellchat <- CellChat::identifyOverExpressedInteractions(oa.cellchat)

oa.cellchat <- CellChat::projectData(oa.cellchat, CellChat::PPI.human)

oa.cellchat <- CellChat::computeCommunProb(oa.cellchat, raw.use = TRUE)

oa.cellchat <- CellChat::filterCommunication(oa.cellchat, min.cells = MIN_CELLS)

oa.cellchat <- CellChat::computeCommunProbPathway(oa.cellchat)

oa.cellchat <- CellChat::aggregateNet(oa.cellchat)

oa.cellchat <- CellChat::netAnalysis_computeCentrality(oa.cellchat, slot.name = "netP")

future::plan("sequential")

group1.net <- CellChat::subsetCommunication(oa.cellchat)

write_csv(group1.net, file.path(TABLE_DIR, "02_group1_net_inter_raw.useT.csv"))

saveRDS(oa.cellchat, file.path(RDS_DIR, "OA_NDUFS3_CellChat_exact_template.rds"))

count.matrix <- oa.cellchat@net$count

weight.matrix <- oa.cellchat@net$weight

write.csv(count.matrix, file.path(TABLE_DIR, "03_interaction_count_matrix.csv"), quote = FALSE)

write.csv(weight.matrix, file.path(TABLE_DIR, "04_interaction_strength_matrix.csv"), quote = FALSE)

