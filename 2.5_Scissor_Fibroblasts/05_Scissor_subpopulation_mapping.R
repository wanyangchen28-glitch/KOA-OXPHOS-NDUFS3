args <- commandArgs(trailingOnly = FALSE)

script_file <- sub("^--file=", "", grep("^--file=", args, value = TRUE)[1])

code_root <- dirname(dirname(normalizePath(script_file, mustWork = TRUE)))

source(file.path(code_root, "00_configuration.R"))

suppressPackageStartupMessages({
    library(Seurat)
    library(dplyr)
})

fib <- readRDS(file.path(input_root, "singlecell", "GSE216651_Fibroblasts_final_7_published_subtypes.rds"))

labels <- read.csv(file.path(output_root, "Scissor", "cell_labels.csv"))

meta <- fib[[]]

meta$cell_id <- rownames(meta)

mapped <- inner_join(meta, labels[, c("cell_id", "Scissor_coefficient", "Scissor_raw_class")], by = "cell_id")

composition <- mapped %>% count(fibroblast_subtype_published, Scissor_raw_class, name = "N")

proportions <- mapped %>% group_by(fibroblast_subtype_published) %>% summarise(N = n(), Scissor_positive = sum(Scissor_raw_class == 
    "Scissor+"), Scissor_positive_proportion = Scissor_positive/N, .groups = "drop")

write.csv(mapped, file.path(output_root, "Scissor", "fibroblast_cell_labels.csv"), row.names = FALSE)

write.csv(composition, file.path(output_root, "Scissor", "fibroblast_composition.csv"), row.names = FALSE)

write.csv(proportions, file.path(output_root, "Scissor", "fibroblast_Scissor_proportions.csv"), row.names = FALSE)

