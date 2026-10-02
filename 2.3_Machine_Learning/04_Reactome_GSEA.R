args <- commandArgs(trailingOnly = FALSE)

script_file <- sub("^--file=", "", grep("^--file=", args, value = TRUE)[1])

code_root <- dirname(dirname(normalizePath(script_file, mustWork = TRUE)))

source(file.path(code_root, "00_configuration.R"))

library(ggplot2)

theme_cns <- theme_bw(base_size = 8, base_family = "Arial") + theme(axis.title = element_text(size = 8), axis.text = element_text(size = 7, 
    color = "#333333"), legend.title = element_text(size = 8, face = "bold"), legend.text = element_text(size = 7), strip.text = element_text(size = 8, 
    face = "bold"), panel.grid = element_blank(), legend.background = element_blank(), legend.key = element_blank())

categorical <- c("#2166AC", "#B2182B", "#1B7837", "#F1A340", "#762A83", "#666666")

categorical_extended <- c("#2166AC", "#B2182B", "#1B7837", "#F1A340", "#762A83", "#666666", "#4393C3", "#D6604D", "#5AAE61", 
    "#B35806", "#9970AB", "#999999")

diverging <- c("#2166AC", "#F7F7F7", "#B2182B")

sequential <- c("#F7FBFF", "#6BAED6", "#08306B")

accent_red <- "#B2182B"

grey <- "#999999"

black <- "#222222"

save_cns_figure <- function(plot, filename, width_mm = 183, height_mm = NULL) {
    ggsave(paste0(filename, ".pdf"), plot, device = cairo_pdf, width = width_mm, height = height_mm, units = "mm", dpi = 300)
    png(paste0(filename, ".png"), width = width_mm, height = height_mm, units = "mm", res = 300, type = "cairo")
    print(plot)
    dev.off()
}

suppressPackageStartupMessages({
    library(clusterProfiler)
    library(dplyr)
    library(ggridges)
    library(readr)
})

project_dir <- file.path(input_root, "bulk_discovery")

expression_file <- file.path(project_dir, "00_preprocess", "06_ComBat_HGNC_maxIQR_expression.csv")

sample_file <- file.path(project_dir, "00_preprocess", "04_sample_info.tsv")

gmt_file <- file.path(project_dir, "08_Hub_gene_GSEA", "00_gene_sets", "ReactomePathways.gmt")

out_dir <- file.path(project_dir, "08_Hub_gene_GSEA", "ATP6V1A_NDUFS3_Reactome_GSEA")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

target_genes <- c("ATP6V1A", "NDUFS3")

seed_value <- 20260916L

min_gs_size <- 10L

max_gs_size <- 500L

top_n_per_gene <- 4L

fdr_threshold <- 0.05

redundancy_threshold <- 0.7

required_files <- c(expression_file, sample_file, gmt_file)

if (any(!file.exists(required_files))) {
    stop("Missing input: ", paste(required_files[!file.exists(required_files)], collapse = "; "))
}

expr_raw <- read.csv(expression_file, check.names = FALSE, stringsAsFactors = FALSE)

gene_symbols <- trimws(as.character(expr_raw[[1]]))

expr <- as.matrix(expr_raw[, -1, drop = FALSE])

storage.mode(expr) <- "double"

rownames(expr) <- gene_symbols

sample_info <- read.delim(sample_file, check.names = FALSE, stringsAsFactors = FALSE)

if (!all(c("GSM", "Group") %in% names(sample_info))) stop("Sample metadata lacks GSM/Group")

common_samples <- intersect(sample_info$GSM, colnames(expr))

expr <- expr[, common_samples, drop = FALSE]

sample_info <- sample_info[match(common_samples, sample_info$GSM), , drop = FALSE]

if (ncol(expr) != 40L) stop("Expected 40 aligned discovery-set samples; found ", ncol(expr))

if (anyDuplicated(rownames(expr))) stop("Expression matrix gene symbols must be unique")

if (!all(target_genes %in% rownames(expr))) stop("Target gene absent from expression matrix")

if (any(!is.finite(expr))) stop("Expression matrix contains non-finite values")

term2gene <- clusterProfiler::read.gmt(gmt_file) %>% transmute(term = term, gene = gene) %>% filter(gene %in% rownames(expr)) %>% 
    distinct(term, gene)

if (n_distinct(term2gene$term) < 100L) stop("Reactome GMT mapping failed")

calculate_ranking <- function(target_gene) {
    target_expression <- as.numeric(expr[target_gene, ])
    rho <- apply(expr, 1, function(x) {
        suppressWarnings(cor(target_expression, x, method = "spearman", use = "pairwise.complete.obs"))
    })
    rho <- rho[is.finite(rho)]
    rho <- rho[names(rho) != target_gene]
    rho_clamped <- pmax(pmin(rho, 0.999999), -0.999999)
    rank_stat <- rho_clamped * sqrt((ncol(expr) - 2)/pmax(1 - rho_clamped^2, 1e-08))
    lexical_order <- order(names(rank_stat))
    tie_break <- numeric(length(rank_stat))
    tie_break[lexical_order] <- seq_along(rank_stat) * 1e-12
    rank_stat <- rank_stat + tie_break
    rank_stat <- sort(rank_stat, decreasing = TRUE)
    correlation_table <- tibble(target_gene = target_gene, gene = names(rank_stat), spearman_rho = unname(rho[names(rank_stat)]), 
        rank_statistic = unname(rank_stat))
    set.seed(seed_value)
    gsea_object <- clusterProfiler::GSEA(geneList = rank_stat, TERM2GENE = term2gene, exponent = 1, minGSSize = min_gs_size, 
        maxGSSize = max_gs_size, eps = 0, pvalueCutoff = 1, pAdjustMethod = "BH", seed = TRUE, by = "fgsea", verbose = FALSE)
    result_table <- as.data.frame(gsea_object@result, stringsAsFactors = FALSE) %>% as_tibble() %>% mutate(target_gene = target_gene, 
        .before = 1)
    list(correlation = correlation_table, result = result_table)
}

analysis_list <- lapply(target_genes, calculate_ranking)

names(analysis_list) <- target_genes

correlation_all <- bind_rows(lapply(analysis_list, `[[`, "correlation"))

gsea_all <- bind_rows(lapply(analysis_list, `[[`, "result"))

if (nrow(gsea_all) == 0L) stop("No Reactome GSEA results were returned")

write_csv(correlation_all, file.path(out_dir, "01_target_gene_Spearman_rankings.csv"))

write_csv(gsea_all, file.path(out_dir, "02_Reactome_GSEA_all_results.csv"))

split_core <- function(x) {
    if (is.na(x) || x == "") 
        character(0)
    else unique(strsplit(x, "/", fixed = TRUE)[[1]])
}

jaccard <- function(a, b) {
    union_n <- length(union(a, b))
    if (union_n == 0L) 
        return(0)
    length(intersect(a, b))/union_n
}

select_nonredundant <- function(df, n_keep = 4L) {
    significant <- df %>% filter(is.finite(NES), is.finite(p.adjust), p.adjust < fdr_threshold, abs(NES) >= 1) %>% arrange(p.adjust, 
        desc(abs(NES)), desc(setSize))
    if (nrow(significant) < n_keep) {
        significant <- df %>% filter(is.finite(NES), is.finite(p.adjust), abs(NES) >= 1) %>% arrange(p.adjust, desc(abs(NES)), 
            desc(setSize))
    }
    chosen <- integer(0)
    chosen_core <- list()
    for (i in seq_len(nrow(significant))) {
        current_core <- split_core(significant$core_enrichment[i])
        overlap <- if (length(chosen_core) == 0L) 
            0
        else max(vapply(chosen_core, function(previous_core) jaccard(current_core, previous_core), numeric(1)))
        if (overlap < redundancy_threshold) {
            chosen <- c(chosen, i)
            chosen_core[[length(chosen_core) + 1L]] <- current_core
        }
        if (length(chosen) >= n_keep) 
            break
    }
    if (length(chosen) < n_keep) {
        extras <- setdiff(seq_len(nrow(significant)), chosen)
        chosen <- c(chosen, head(extras, n_keep - length(chosen)))
    }
    significant[chosen, , drop = FALSE]
}

selected_stats <- gsea_all %>% group_by(target_gene) %>% group_split() %>% lapply(select_nonredundant, n_keep = top_n_per_gene) %>% 
    bind_rows() %>% mutate(FDR = p.adjust, neglog10_fdr = pmin(-log10(pmax(FDR, .Machine$double.xmin)), 12), abs_NES = abs(NES))

if (any(table(selected_stats$target_gene) != top_n_per_gene)) {
    stop("Could not select four pathways for each target gene")
}

wrap_label <- function(x, width = 36L) {
    vapply(x, function(z) paste(strwrap(z, width = width), collapse = "\n"), character(1))
}

selected_stats <- selected_stats %>% group_by(target_gene) %>% arrange(desc(NES), .by_group = TRUE) %>% ungroup() %>% mutate(compact_description = recode(Description, 
    `Aerobic respiration and respiratory electron transport` = "Aerobic respiration & electron transport", `Nucleotide Excision Repair` = "Nucleotide excision repair", 
    .default = Description), pathway_label = wrap_label(compact_description), display_order = row_number())

write_csv(selected_stats, file.path(out_dir, "03_selected_pathway_statistics.csv"))

ridge_df <- selected_stats %>% select(target_gene, Description, pathway_label, NES, FDR, neglog10_fdr, abs_NES) %>% inner_join(term2gene, 
    by = c(Description = "term"), relationship = "many-to-many") %>% inner_join(correlation_all, by = c("target_gene", "gene"))

selected_counts <- ridge_df %>% count(target_gene, Description, name = "mapped_gene_count")

selected_stats <- selected_stats %>% left_join(selected_counts, by = c("target_gene", "Description"))

write_csv(selected_stats, file.path(out_dir, "03_selected_pathway_statistics.csv"))

write_csv(ridge_df, file.path(out_dir, "04_selected_pathway_member_rank_statistics.csv"))

if (any(selected_stats$mapped_gene_count < min_gs_size)) stop("Selected pathway has too few mapped genes")

level_order <- rev(unique(selected_stats$pathway_label))

ridge_df <- ridge_df %>% mutate(pathway_label = factor(pathway_label, levels = level_order))

selected_stats <- selected_stats %>% mutate(pathway_label = factor(pathway_label, levels = level_order))

rank_quantiles <- quantile(ridge_df$rank_statistic, probs = c(0.005, 0.995), na.rm = TRUE)

rank_span <- diff(rank_quantiles)

x_dot <- unname(rank_quantiles[1] - 0.16 * rank_span)

x_limits <- c(x_dot - 0.08 * rank_span, unname(rank_quantiles[2] + 0.05 * rank_span))

p <- ggplot(ridge_df, aes(x = rank_statistic, y = pathway_label, group = interaction(target_gene, Description), fill = neglog10_fdr)) + 
    geom_vline(xintercept = 0, color = "#666666", linewidth = 0.44) + geom_density_ridges(scale = 0.78, rel_min_height = 0.001, 
    bandwidth = 0.42, color = "#666666", linewidth = 0.38, alpha = 0.94, na.rm = TRUE) + geom_point(aes(y = pathway_label), 
    color = "#3F3F3F", shape = 124, size = 1.45, alpha = 0.38, show.legend = FALSE) + geom_point(data = selected_stats, aes(x = x_dot, 
    y = pathway_label, color = NES), inherit.aes = FALSE, shape = 16, size = 4.2) + scale_fill_gradientn(colours = c("#67A9CF", 
    "#A6DBA0", "#FEE08B", "#F46D43", "#D73027"), name = expression(-log[10](FDR)), guide = guide_colorbar(order = 2, direction = "vertical", 
    title.position = "top", barwidth = grid::unit(3.5, "mm"), barheight = grid::unit(28, "mm"))) + scale_color_gradient2(low = "#2C7BB6", 
    mid = "#F7F7F7", high = "#D7191C", midpoint = 0, name = "NES", guide = guide_colorbar(order = 1, direction = "vertical", 
        title.position = "top", barwidth = grid::unit(3.5, "mm"), barheight = grid::unit(28, "mm"))) + scale_y_discrete(expand = expansion(add = c(0.16, 
    0.82))) + facet_grid(rows = vars(target_gene), scales = "free_y", space = "free_y", switch = "y") + coord_cartesian(xlim = x_limits, 
    clip = "off") + labs(title = "Reactome GSEA ridge plot", subtitle = NULL, x = "GSEA rank metric", y = NULL, caption = NULL) + 
    theme_cns + theme(text = element_text(family = "Helvetica"), panel.background = element_rect(fill = "white", color = NA), 
    panel.border = element_blank(), panel.grid = element_blank(), axis.line.x = element_line(color = black, linewidth = 0.52), 
    axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_text(size = 12.5, color = black, 
        lineheight = 0.92, hjust = 1), axis.text.x = element_text(size = 10.5, color = black), axis.title.x = element_text(size = 12, 
        margin = margin(t = 5)), plot.title = element_text(size = 15, face = "bold", hjust = 0, margin = margin(b = 5)), 
    strip.background = element_blank(), strip.placement = "outside", strip.text.y.left = element_text(size = 12.5, face = "bold", 
        angle = 90, color = black), panel.spacing.y = grid::unit(0.6, "mm"), legend.position = "right", legend.box = "vertical", 
    legend.direction = "vertical", legend.title = element_text(size = 11.5, face = "bold"), legend.text = element_text(size = 10.5), 
    legend.spacing.y = grid::unit(2.5, "mm"), legend.margin = margin(l = 2, unit = "mm"), plot.margin = margin(3, 3, 3, 3, 
        unit = "mm"))

output_stem <- file.path(out_dir, "Fig2_ATP6V1A_NDUFS3_Reactome_GSEA_ridgeplot")

save_cns_figure(p, output_stem, width_mm = 183, height_mm = 145)

invisible(NULL)

writeLines(capture.output(sessionInfo()), file.path(out_dir, "06_sessionInfo.txt"))

message("Completed: ", output_stem, ".pdf and .png")

