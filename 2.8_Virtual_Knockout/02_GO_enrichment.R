args <- commandArgs(trailingOnly = FALSE)

script_file <- sub("^--file=", "", grep("^--file=", args, value = TRUE)[1])

code_root <- dirname(dirname(normalizePath(script_file, mustWork = TRUE)))

source(file.path(code_root, "00_configuration.R"))

gc()

suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(ragg)
    library(Cairo)
    library(clusterProfiler)
    library(org.Hs.eg.db)
    library(AnnotationDbi)
})

diff_file <- paste0(paste0(file.path(input_root, "external_assets/PRG4_lining_virtual_knockout"), "/"), "OA_PRG4_lining_fibroblasts_NDUFS3_KO_diffRegulation_all_genes.csv")

output_dir <- paste0(paste0(file.path(output_root, "outputs"), "/"), "45_PRG4_NDUFS3_KO_GO_exact_Figure7G_style")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

font_family <- "Arial"

target_gene <- "NDUFS3"

theme_large <- function(base_size = 16) {
    theme_classic(base_size = base_size, base_family = font_family) + theme(text = element_text(family = font_family, colour = "black"), 
        axis.text = element_text(family = font_family, colour = "black", size = base_size - 1), axis.title = element_text(family = font_family, 
            colour = "black", size = base_size + 0.5), plot.title = element_text(family = font_family, face = "bold", size = base_size + 
            2, hjust = 0.5, margin = margin(b = 5)), plot.subtitle = element_text(family = font_family, size = base_size - 
            2.5, hjust = 0.5, colour = "grey30", margin = margin(b = 6)), legend.title = element_text(family = font_family, 
            face = "bold", size = base_size - 2), legend.text = element_text(family = font_family, size = base_size - 2.5), 
        axis.line = element_line(linewidth = 0.55, colour = "black"), axis.ticks = element_line(linewidth = 0.45, colour = "black"), 
        axis.ticks.length = grid::unit(2.4, "mm"))
}

save_panel <- function(plot_object, stem, width = 11, height = 7.8, dpi = 500) {
    pdf_file <- file.path(output_dir, paste0(stem, ".pdf"))
    png_file <- file.path(output_dir, paste0(stem, ".png"))
    svg_file <- file.path(output_dir, paste0(stem, ".svg"))
    while (grDevices::dev.cur() != 1L) grDevices::dev.off()
    Cairo::CairoPDF(file = pdf_file, width = width, height = height, family = font_family, onefile = TRUE, bg = "white")
    print(plot_object)
    grDevices::dev.off()
    ragg::agg_png(filename = png_file, width = width, height = height, units = "in", res = dpi, background = "white", scaling = 1)
    print(plot_object)
    grDevices::dev.off()
    svglite::svglite(svg_file, width = width, height = height)
    print(plot_object)
    grDevices::dev.off()
}

diff_df <- read.csv(diff_file, stringsAsFactors = FALSE, check.names = FALSE)

genes_use <- unique(pull(filter(diff_df, is.finite(as.numeric(p.adj)), as.numeric(p.adj) < 0.05, gene != target_gene), gene))

gene_map <- distinct(filter(AnnotationDbi::select(org.Hs.eg.db, keys = genes_use, columns = c("SYMBOL", "ENTREZID"), keytype = "SYMBOL"), 
    !is.na(ENTREZID)), SYMBOL, ENTREZID)

run_go <- function(ontology) {
    ego <- clusterProfiler::enrichGO(gene = unique(gene_map$ENTREZID), OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = ontology, 
        pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = TRUE)
    mutate(as.data.frame(ego), ONTOLOGY = ontology, Count = as.numeric(Count), p.adjust = as.numeric(p.adjust), neg_log10_fdr = -log10(p.adjust))
}

go_bp <- run_go("BP")

go_cc <- run_go("CC")

go_mf <- run_go("MF")

write.csv(go_bp, file.path(output_dir, "GO_BP.csv"), row.names = FALSE)

write.csv(go_cc, file.path(output_dir, "GO_CC.csv"), row.names = FALSE)

write.csv(go_mf, file.path(output_dir, "GO_MF.csv"), row.names = FALSE)

