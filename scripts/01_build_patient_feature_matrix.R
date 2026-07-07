############################################################
# 01_build_patient_feature_matrix.R
#
# Project:
# AI-Driven Patient Stratification and Precision Targeting
# in Lupus Nephritis
#
# Purpose:
# Build a patient-level feature matrix from the AMP lupus
# nephritis single-cell atlas.
############################################################

library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(tibble)
library(pheatmap)

############################################################
# Load data
############################################################

ln <- readRDS(
  "data/AMP_LN_human_final.rds"
)

meta <- read.delim(
  "data/celseq_meta.tsv"
)

############################################################
# Add sample metadata to Seurat object
############################################################

meta_small <- meta %>%
  select(
    cell_name,
    sample,
    disease,
    type
  )

rownames(meta_small) <- meta_small$cell_name

ln <- AddMetaData(
  ln,
  metadata = meta_small
)

ln$celltype <- Idents(ln)

############################################################
# Patient-level cell composition
############################################################

cell_composition <- ln@meta.data %>%
  filter(
    sample != "none",
    !grepl("Double", sample)
  ) %>%
  count(
    sample,
    disease,
    celltype
  ) %>%
  group_by(
    sample,
    disease
  ) %>%
  mutate(
    fraction = n / sum(n)
  ) %>%
  ungroup() %>%
  select(
    sample,
    disease,
    celltype,
    fraction
  ) %>%
  pivot_wider(
    names_from = celltype,
    values_from = fraction,
    values_fill = 0
  )

write.csv(
  cell_composition,
  "results/Patient_Cell_Composition.csv",
  row.names = FALSE
)

############################################################
# Patient-level target expression
############################################################

target_genes <- c(
  "C5AR1",
  "CSF1R",
  "LILRB2",
  "PILRA",
  "SIGLEC1",
  "CLEC7A",
  "TLR4",
  "P2RX7",
  "C3AR1",
  "CD300E"
)

target_expr <- AverageExpression(
  ln,
  features = target_genes,
  group.by = "sample"
)$RNA

target_expr <- as.data.frame(
  t(target_expr)
)

target_expr$sample <- rownames(target_expr)

target_expr <- target_expr %>%
  filter(
    sample != "none",
    !grepl("Double", sample)
  )

target_expr$sample <- sub(
  "^g",
  "",
  target_expr$sample
)

############################################################
# Build patient feature matrix
############################################################

patient_features <- cell_composition %>%
  left_join(
    target_expr,
    by = "sample"
  )

write.csv(
  patient_features,
  "results/Patient_Feature_Matrix_v1.csv",
  row.names = FALSE
)

############################################################
# Patient feature heatmap
############################################################

patient_matrix <- patient_features %>%
  column_to_rownames("sample") %>%
  select(-disease)

png(
  "figures/Figure1_Patient_Feature_Heatmap.png",
  width = 2600,
  height = 2200,
  res = 300
)

pheatmap(
  patient_matrix,
  scale = "column",
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "ward.D2",
  main = "Patient-level single-cell and target feature matrix"
)

dev.off()

############################################################
# Output checks
############################################################

dim(patient_features)
list.files("figures")
list.files("results")