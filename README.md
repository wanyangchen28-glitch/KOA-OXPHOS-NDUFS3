# KOA-OXPHOS-NDUFS3

Integrative Analysis of Bulk and Single-Cell Transcriptomes with Machine Learning Identifies NDUFS3 as a Candidate Gene Associated with OXPHOS Dysregulation in Knee Osteoarthritis Synovium.

## Folder Contents

- `2.1_Data_Sources`: Data acquisition and preprocessing.
- `2.2_DEGs_WGCNA_GSEA`: Differential expression analysis, WGCNA, and GSEA.
- `2.3_Machine_Learning`: Feature selection, model development, and external validation.
- `2.4_Single_Cell_Atlas`: Single-cell quality control, integration, clustering, and annotation.
- `2.5_Scissor_Fibroblasts`: Scissor analysis and fibroblast subpopulation analysis.
- `2.6_OXPHOS_Scoring`: OXPHOS scoring and candidate gene expression analysis.
- `2.7_CellChat_Pseudotime`: Cell–cell communication and pseudotime analysis.
- `2.8_Virtual_Knockout`: NDUFS3 virtual knockout analysis.
-  `2.9_Drug_Screening`: DrugReflector-based candidate compound screening, molecular docking, and molecular dynamics simulations.
- `data`: Sample metadata.
- `00_configuration.R`: Shared configuration.

## Datasets

GEO accession numbers: **GSE55235, GSE55457, GSE32317, GSE89408, and GSE216651**.

## Software Requirements

R (version 4.5.2) and Python.

Main R packages include limma, sva, WGCNA, fgsea, glmnet, randomForest, Boruta, caret, e1071, pROC, clusterProfiler, Seurat, scDblFinder, Scissor, AUCell, UCell, singscore, GSVA, CellChat, Monocle 2, scTenifoldKnk, ggplot2, and ComplexHeatmap. DrugReflector is used for candidate compound screening.
Molecular docking and molecular dynamics analyses require fpocket, AutoDock Vina, PyMOL, AmberTools, and GROMACS (version 2025.4).
