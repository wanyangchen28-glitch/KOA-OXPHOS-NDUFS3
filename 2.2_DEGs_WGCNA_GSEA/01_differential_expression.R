args <- commandArgs(trailingOnly = FALSE)

script_file <- sub("^--file=", "", grep("^--file=", args, value = TRUE)[1])

code_root <- dirname(dirname(normalizePath(script_file, mustWork = TRUE)))

source(file.path(code_root, "00_configuration.R"))

gc()

options(digits = 4)

base_dir <- normalizePath(file.path(input_root, "bulk_discovery"), mustWork = TRUE)

work_dir <- file.path(base_dir, "01_DEGs")

if (!dir.exists(work_dir)) {
    dir.create(work_dir, recursive = TRUE)
}

setwd(work_dir)

logFC_threshold <- 0.5

p_adj_threshold <- 0.05

lfc_line <- logFC_threshold

padj_line <- p_adj_threshold

ctrl_name <- "Control"

case_name <- "OA"

group_levels <- c(ctrl_name, case_name)

suppressPackageStartupMessages({
    library(magrittr)
    library(limma)
    library(dplyr)
    library(data.table)
    library(ggplot2)
    library(ggthemes)
    library(ggrepel)
    library(ComplexHeatmap)
    library(circlize)
    library(grid)
})

genes_expr <- read.csv(file.path(base_dir, "expression_matrix.csv"), header = TRUE, row.names = 1, check.names = FALSE)

group_list <- read.csv(file.path(base_dir, "group_metadata.csv"), header = TRUE, stringsAsFactors = FALSE)

if (!all(c("Sample", "Group") %in% colnames(group_list))) {
    stop("group_metadata.csv 必须包含 Sample 和 Group 两列。")
}

rownames(group_list) <- group_list$Sample

group_list$Group <- factor(group_list$Group, levels = group_levels)

if (any(is.na(group_list$Group))) {
    stop("group_metadata.csv 中存在不属于 ctrl_name 或 case_name 的分组，请检查 Group 列。")
}

common_samples <- intersect(colnames(genes_expr), rownames(group_list))

if (length(common_samples) == 0) {
    stop("表达矩阵列名和 metadata 的 Sample 没有匹配样本。")
}

genes_expr <- genes_expr[, common_samples, drop = FALSE]

group_list <- group_list[common_samples, , drop = FALSE]

stopifnot(all(colnames(genes_expr) == rownames(group_list)))

genes_expr <- as.matrix(genes_expr)

mode(genes_expr) <- "numeric"

design <- model.matrix(~0 + Group, data = group_list)

design_names <- make.names(levels(group_list$Group))

colnames(design) <- design_names

rownames(design) <- colnames(genes_expr)

ctrl_design <- make.names(ctrl_name)

case_design <- make.names(case_name)

contrast_str <- paste0(case_design, "-", ctrl_design)

contrast.matrix <- makeContrasts(contrasts = contrast_str, levels = design)

fit <- lmFit(genes_expr, design)

fit2 <- contrasts.fit(fit, contrast.matrix)

fit2 <- eBayes(fit2)

DEG <- topTable(fit2, coef = 1, n = Inf, adjust.method = "BH")

colnames(DEG)[colnames(DEG) == "logFC"] <- "log2FoldChange"

DEG$change <- ifelse(DEG$adj.P.Val >= p_adj_threshold, "Not", ifelse(DEG$log2FoldChange >= logFC_threshold, "Up", ifelse(DEG$log2FoldChange <= 
    -logFC_threshold, "Down", "Not")))

DEG$gene <- rownames(DEG)

DEG$Symbols <- rownames(DEG)

write.csv(DEG, file = "KOA_GPL96_merged_all_DEGs_fc.csv", quote = FALSE, row.names = TRUE)

DEG_sig <- DEG %>% dplyr::filter(adj.P.Val < p_adj_threshold, abs(log2FoldChange) >= logFC_threshold)

write.csv(DEG_sig, file = "KOA_GPL96_merged_sig_DEGs_fc.csv", quote = FALSE, row.names = TRUE)

data.table::fwrite(as.data.frame(DEG_sig[DEG_sig$change == "Down", ]), file = "KOA_GPL96_merged_DEG_FC_Down.csv", row.names = TRUE)

data.table::fwrite(as.data.frame(DEG_sig[DEG_sig$change == "Up", ]), file = "KOA_GPL96_merged_DEG_FC_Up.csv", row.names = TRUE)

positive_p <- DEG$adj.P.Val[DEG$adj.P.Val > 0]

if (length(positive_p) == 0) {
    stop("adj.P.Val 全部为 0 或缺失，无法绘制火山图。")
}

min_positive_p <- min(positive_p, na.rm = TRUE)

volcano_data <- DEG %>% dplyr::mutate(gene_symbol = gene, logFC = log2FoldChange, plot_padj = adj.P.Val, plot_padj = ifelse(is.na(plot_padj), 
    NA, plot_padj), plot_padj = ifelse(plot_padj <= 0, min_positive_p, plot_padj), neg_log10_padj = -log10(plot_padj), change = dplyr::case_when(logFC >= 
    lfc_line & plot_padj < padj_line ~ "Up", logFC <= -lfc_line & plot_padj < padj_line ~ "Down", TRUE ~ "Not")) %>% dplyr::filter(!is.na(logFC), 
    !is.na(plot_padj), is.finite(neg_log10_padj))

volcano_up <- volcano_data %>% dplyr::filter(change == "Up")

volcano_down <- volcano_data %>% dplyr::filter(change == "Down")

volcano_sig <- volcano_data %>% dplyr::filter(change != "Not")

n_up <- nrow(volcano_up)

n_down <- nrow(volcano_down)

n_sig <- nrow(volcano_sig)

label_up <- volcano_up %>% dplyr::arrange(plot_padj, dplyr::desc(abs(logFC))) %>% dplyr::slice_head(n = 10)

label_down <- volcano_down %>% dplyr::arrange(plot_padj, dplyr::desc(abs(logFC))) %>% dplyr::slice_head(n = 10)

data_repel <- rbind(label_up, label_down)

control_group <- ctrl_name

case_group <- case_name

if (requireNamespace("ggtext", quietly = TRUE)) {
    plot_title <- paste0("<span style='color:#d73027;'>", case_group, "</span>", " vs ", "<span style='color:#3288bd;'>", 
        control_group, "</span>")
    title_theme <- ggtext::element_markdown(hjust = 0.5, face = "bold", size = 22, margin = margin(b = 6))
} else {
    plot_title <- paste(case_group, "vs", control_group)
    title_theme <- element_text(hjust = 0.5, face = "bold", size = 22, margin = margin(b = 6))
}

plot_subtitle <- paste0("Adj.P.Value = ", padj_line, "; ", "log2FC = ", lfc_line, "; ", "Up: ", n_up, "; ", "Down: ", n_down, 
    "; ", "Total: ", n_sig)

y_max <- max(volcano_data$neg_log10_padj, na.rm = TRUE)

y_lim <- ceiling(y_max + 2)

x_abs <- max(abs(volcano_data$logFC), na.rm = TRUE)

x_lim <- max(0.3, ceiling(x_abs * 10)/10)

x_min <- -x_lim

x_max <- x_lim

arrow_y <- y_lim * 0.94

volcano_plot2 <- ggplot(volcano_data, aes(x = logFC, y = neg_log10_padj)) + geom_point(aes(color = logFC, size = neg_log10_padj), 
    alpha = 0.8, na.rm = TRUE) + geom_text_repel(data = label_up, aes(label = gene_symbol), color = "#d73027", segment.color = "#d73027", 
    box.padding = 0.5, nudge_x = 0.15, nudge_y = 0.2, segment.curvature = -0.1, segment.ncp = 3, segment.angle = 10, direction = "y", 
    hjust = "left", max.overlaps = 200, family = "Times", size = 4, show.legend = FALSE) + geom_text_repel(data = label_down, 
    aes(label = gene_symbol), color = "#3288bd", segment.color = "#3288bd", box.padding = 0.5, nudge_x = -0.15, nudge_y = 0.2, 
    segment.curvature = -0.1, segment.ncp = 3, segment.angle = 20, direction = "y", hjust = "right", max.overlaps = 200, 
    family = "Times", size = 4, show.legend = FALSE) + scale_color_gradientn(colours = c("#3288bd", "#66c2a5", "#ffffbf", 
    "#f46d43", "#9e0142"), name = "log2FC") + scale_size(range = c(1, 7), name = "-log10(Adj.P.Value)") + geom_vline(xintercept = c(-lfc_line, 
    lfc_line), linetype = 4, color = "darkgray", linewidth = 0.6) + geom_hline(yintercept = -log10(padj_line), linetype = 4, 
    color = "darkgray", linewidth = 0.6) + annotate("segment", x = -lfc_line, xend = x_min + 0.15, y = arrow_y, yend = arrow_y, 
    arrow = arrow(angle = 45, length = unit(0.2, "cm"), ends = "last"), linewidth = 1, color = "#74add1") + annotate("text", 
    x = (x_min + (-lfc_line))/2, y = arrow_y, label = "DOWN", color = "#74add1", vjust = -0.8, size = 4, family = "Times", 
    fontface = "bold") + annotate("segment", x = lfc_line, xend = x_max - 0.15, y = arrow_y, yend = arrow_y, arrow = arrow(angle = 45, 
    length = unit(0.2, "cm"), ends = "last"), linewidth = 1, color = "#d73027") + annotate("text", x = (x_max + lfc_line)/2, 
    y = arrow_y, label = "UP", color = "#d73027", vjust = -0.8, size = 4, family = "Times", fontface = "bold") + labs(x = "log2 (Fold Change)", 
    y = "-log10 (Adj.P.Value)", title = plot_title, subtitle = plot_subtitle) + coord_cartesian(xlim = c(x_min, x_max), ylim = c(0, 
    y_lim), clip = "off") + theme_bw(base_family = "Times") + theme(panel.grid = element_blank(), legend.background = element_rect(color = "#808080", 
    linetype = 1), legend.title = element_text(face = "bold", size = 11), legend.text = element_text(size = 10), axis.text = element_text(size = 13, 
    color = "#000000", face = "bold"), axis.title = element_text(size = 15, face = "bold"), plot.title = title_theme, plot.subtitle = element_text(hjust = 0.5, 
    size = 12, color = "black", face = "italic", margin = margin(b = 10)), plot.margin = margin(t = 20, r = 40, b = 10, l = 20))

ggsave(filename = "Volcano_plot.png", plot = volcano_plot2, width = 10, height = 8, dpi = 600)

ggsave(filename = "Volcano_plot.pdf", plot = volcano_plot2, width = 10, height = 8, family = "Times")

if (nrow(data_repel) > 0) {
    heat_genes <- unique(data_repel$gene_symbol)
} else {
    heat_genes <- DEG %>% dplyr::arrange(adj.P.Val, dplyr::desc(abs(log2FoldChange))) %>% dplyr::slice_head(n = 20) %>% dplyr::pull(gene)
}

heat_genes <- intersect(heat_genes, rownames(genes_expr))

if (length(heat_genes) < 2) {
    stop("可用于热图的基因少于 2 个，请检查差异分析结果或基因名是否匹配。")
}

mat <- genes_expr[heat_genes, , drop = FALSE]

gene_order <- DEG[heat_genes, , drop = FALSE] %>% dplyr::arrange(dplyr::desc(log2FoldChange)) %>% rownames()

mat <- mat[gene_order, , drop = FALSE]

mat_scaled <- t(scale(t(mat)))

mat_scaled[is.na(mat_scaled)] <- 0

mat_scaled[mat_scaled < -2] <- -2

mat_scaled[mat_scaled > 2] <- 2

heatmap_colors <- list(Group = setNames(c("#2187f4", "#ff3030"), group_levels))

draw_density_heatmap <- function() {
    densityHeatmap(mat_scaled, title = "Distribution as heatmap", ylab = " ", height = unit(6, "cm")) %v% HeatmapAnnotation(Group = factor(group_list$Group, 
        levels = group_levels), col = heatmap_colors) %v% Heatmap(mat_scaled, row_names_gp = gpar(fontsize = 9, fontfamily = "Times"), 
        column_names_gp = gpar(fontsize = 9, fontfamily = "Times"), show_column_names = FALSE, show_row_names = TRUE, name = "expression", 
        cluster_rows = FALSE, cluster_columns = FALSE, height = unit(6, "cm"), col = circlize::colorRamp2(c(-2, 0, 2), c("deepskyblue", 
            "white", "deeppink")))
}

png(filename = "Density_heatmap.png", width = 10, height = 8, units = "in", res = 600, family = "Times")

draw_density_heatmap()

dev.off()

pdf(file = "Density_heatmap.pdf", width = 10, height = 8, family = "Times")

draw_density_heatmap()

dev.off()

cat("Analysis Completed!\n")

