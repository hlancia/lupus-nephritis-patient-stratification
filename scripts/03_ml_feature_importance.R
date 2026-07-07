############################################################
# 03_ml_feature_importance.R
#
# Project:
# AI-Driven Patient Stratification and Precision Targeting
# in Lupus Nephritis
#
# Purpose:
# Use Random Forest-style exploratory machine learning to
# identify features that drive patient cluster assignment.
############################################################

library(dplyr)
library(ggplot2)
library(tibble)
library(randomForest)

############################################################
# Load patient stratification table
############################################################

patient_data <- read.csv(
  "results/Patient_Stratification_Table.csv"
)

############################################################
# Prepare ML matrix
############################################################

ml_data <- patient_data %>%
  select(
    -sample,
    -disease
  )

ml_data$Patient_Cluster <- as.factor(
  ml_data$Patient_Cluster
)


############################################################
# Train exploratory Random Forest model
############################################################

set.seed(123)

rf_model <- randomForest(
  Patient_Cluster ~ .,
  data = ml_data,
  importance = TRUE,
  ntree = 1000
)

rf_model


############################################################
# Extract feature importance
############################################################

importance_df <- as.data.frame(
  importance(rf_model)
)

importance_df$Feature <- rownames(importance_df)

importance_df <- importance_df %>%
  arrange(desc(MeanDecreaseGini))

write.csv(
  importance_df,
  "results/ML_Feature_Importance.csv",
  row.names = FALSE
)


############################################################
# Plot feature importance
############################################################

png(
  "figures/Figure5_ML_Feature_Importance.png",
  width = 2200,
  height = 1600,
  res = 300
)

ggplot(
  head(importance_df, 15),
  aes(
    x = reorder(Feature, MeanDecreaseGini),
    y = MeanDecreaseGini
  )
) +
  geom_col() +
  coord_flip() +
  theme_classic() +
  labs(
    title = "Features driving patient subtype assignment",
    x = "Feature",
    y = "Mean decrease in Gini"
  )

dev.off()

############################################################
# Output checks
############################################################

head(importance_df, 20)

list.files("figures")
list.files("results")