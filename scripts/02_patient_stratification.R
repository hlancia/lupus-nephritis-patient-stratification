############################################################
# 02_patient_stratification.R
#
# Project:
# AI-Driven Patient Stratification and Precision Targeting
# in Lupus Nephritis
#
# Purpose:
# Perform unsupervised patient stratification using the
# patient-level single-cell and target feature matrix.
############################################################

library(dplyr)
library(ggplot2)
library(tibble)
library(pheatmap)

############################################################
# Load patient feature matrix
############################################################

patient_features <- read.csv(
  "results/Patient_Feature_Matrix_v1.csv"
)

patient_matrix <- patient_features %>%
  column_to_rownames("sample") %>%
  select(-disease)

patient_scaled <- scale(patient_matrix)


############################################################
# PCA
############################################################

pca_res <- prcomp(
  patient_scaled,
  center = TRUE,
  scale. = FALSE
)

pca_df <- data.frame(
  sample = rownames(patient_scaled),
  disease = patient_features$disease,
  PC1 = pca_res$x[, 1],
  PC2 = pca_res$x[, 2]
)

png(
  "figures/Figure2_Patient_PCA.png",
  width = 2200,
  height = 1600,
  res = 300
)

ggplot(
  pca_df,
  aes(
    x = PC1,
    y = PC2,
    color = disease,
    label = sample
  )
) +
  geom_point(size = 4) +
  geom_text(vjust = -0.8, size = 3) +
  theme_classic() +
  labs(
    title = "PCA-based patient stratification in lupus nephritis",
    x = "PC1",
    y = "PC2"
  )

dev.off()


############################################################
# Hierarchical clustering
############################################################

patient_dist <- dist(patient_scaled)

patient_hc <- hclust(
  patient_dist,
  method = "ward.D2"
)

png(
  "figures/Figure3_Patient_Hierarchical_Clustering.png",
  width = 2400,
  height = 1600,
  res = 300
)

plot(
  patient_hc,
  main = "Hierarchical clustering of lupus nephritis patients",
  xlab = "",
  sub = ""
)

dev.off()


############################################################
# Define patient clusters
############################################################

patient_clusters <- cutree(
  patient_hc,
  k = 3
)

cluster_table <- data.frame(
  sample = names(patient_clusters),
  Patient_Cluster = paste0("Cluster_", patient_clusters)
)

patient_stratification <- patient_features %>%
  left_join(
    cluster_table,
    by = "sample"
  )

write.csv(
  patient_stratification,
  "results/Patient_Stratification_Table.csv",
  row.names = FALSE
)

write.csv(
  pca_df,
  "results/Patient_PCA_Coordinates.csv",
  row.names = FALSE
)


############################################################
# Heatmap with patient clusters
############################################################

annotation_row <- data.frame(
  Cluster = patient_stratification$Patient_Cluster,
  Disease = patient_stratification$disease
)

rownames(annotation_row) <- patient_stratification$sample

png(
  "figures/Figure4_Patient_Cluster_Heatmap.png",
  width = 2600,
  height = 2200,
  res = 300
)

pheatmap(
  patient_scaled,
  annotation_row = annotation_row,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "ward.D2",
  main = "Patient clusters based on cell composition and target activity"
)

dev.off()


############################################################
# Output checks
############################################################

table(patient_stratification$Patient_Cluster)

list.files("figures")
list.files("results")