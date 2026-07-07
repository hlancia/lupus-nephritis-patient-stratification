############################################################
# 04_therapeutic_recommendation_engine.R
#
# Project:
# AI-Driven Patient Stratification and Precision Targeting
# in Lupus Nephritis
#
# Purpose:
# Map patient clusters to prioritized target programs and
# generate interpretable precision therapeutic hypotheses.
############################################################

library(dplyr)
library(tidyr)
library(ggplot2)
library(pheatmap)

############################################################
# Load data
############################################################

patient_data <- read.csv(
  "results/Patient_Stratification_Table.csv"
)

feature_importance <- read.csv(
  "results/ML_Feature_Importance.csv"
)

############################################################
# Define therapeutic target genes
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

############################################################
# Cluster-level target activity
############################################################

cluster_target_activity <- patient_data %>%
  group_by(Patient_Cluster) %>%
  summarise(
    across(
      all_of(target_genes),
      mean,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

write.csv(
  cluster_target_activity,
  "results/Cluster_Target_Activity.csv",
  row.names = FALSE
)

############################################################
# Long-format cluster target ranking
############################################################

cluster_target_ranking <- cluster_target_activity %>%
  pivot_longer(
    cols = all_of(target_genes),
    names_to = "Target",
    values_to = "Mean_Activity"
  ) %>%
  group_by(Patient_Cluster) %>%
  arrange(desc(Mean_Activity), .by_group = TRUE) %>%
  mutate(
    Cluster_Target_Rank = row_number()
  ) %>%
  ungroup()

write.csv(
  cluster_target_ranking,
  "results/Cluster_Target_Ranking.csv",
  row.names = FALSE
)

############################################################
# Heatmap of target activity by cluster
############################################################

cluster_matrix <- cluster_target_activity
rownames(cluster_matrix) <- cluster_matrix$Patient_Cluster
cluster_matrix$Patient_Cluster <- NULL

png(
  "figures/Figure6_Cluster_Target_Programs.png",
  width = 2200,
  height = 1600,
  res = 300
)

pheatmap(
  cluster_matrix,
  scale = "column",
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "ward.D2",
  main = "Cluster-specific therapeutic target programs"
)

dev.off()

############################################################
# Therapeutic hypothesis table
############################################################

top_targets_per_cluster <- cluster_target_ranking %>%
  filter(Cluster_Target_Rank <= 3)

therapeutic_hypotheses <- top_targets_per_cluster %>%
  mutate(
    Therapeutic_Hypothesis = case_when(
      Target == "C5AR1" ~ "Complement/C5aR1-high myeloid inflammatory program",
      Target == "CSF1R" ~ "Macrophage survival and differentiation program",
      Target == "LILRB2" ~ "Inhibitory myeloid immune checkpoint program",
      Target == "PILRA" ~ "Myeloid immune regulatory receptor program",
      Target == "SIGLEC1" ~ "Interferon-associated macrophage activation program",
      Target == "CLEC7A" ~ "Pattern-recognition and innate immune activation program",
      Target == "TLR4" ~ "TLR-driven inflammatory activation program",
      Target == "P2RX7" ~ "Inflammasome/purinergic inflammatory signaling program",
      Target == "C3AR1" ~ "Complement/C3aR1-associated myeloid program",
      Target == "CD300E" ~ "Myeloid activation receptor program",
      TRUE ~ "Target-associated immune program"
    )
  )

write.csv(
  therapeutic_hypotheses,
  "results/Therapeutic_Hypotheses_By_Cluster.csv",
  row.names = FALSE
)

############################################################
# Output checks
############################################################

cluster_target_activity

top_targets_per_cluster

therapeutic_hypotheses

list.files("figures")
list.files("results")