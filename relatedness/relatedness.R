## Script to process and analyze relatedness data
#   relatedness was calculated on Beartooth using vcftools

# load required packages
library(tidyverse)
library(vegan)
library(Matrix)

# setwd("C:/Users/UW-User/Desktop/")
setwd("~/UWyo/kellyn_chipmunks/Analyses_git/relatedness/") # Sean's working dir

# Read the relatedness data from the file
relatedness_data <- read.table("chip_c93_branch1.relatedness", header = TRUE)

# Create a matrix using pivot_wider
relatedness_matrix <- relatedness_data %>%
  pivot_wider(names_from = INDV2, values_from = RELATEDNESS_AJK) %>%
  column_to_rownames(var = "INDV1")

# Print the resulting matrix
print(relatedness_matrix)



### Weight data

# Read the Euthanasia Log data from the file
euthanasia_data <- read.csv("Euthanasia_Log.csv")

# Remove duplicates based on "Chipmunk.ID"
euthanasia_data <- euthanasia_data %>% distinct(Chipmunk.ID, .keep_all = TRUE)


# make weight column numeric
euthanasia_data$Weight..g. <- as.numeric(euthanasia_data$Weight..g.)

# Remove NAs
euthanasia_data <- euthanasia_data[!is.na(euthanasia_data$Weight..g.),]

# Extract Chipmunk.ID and Weight..g. columns
euthanasia_subset <- euthanasia_data[, c("Chipmunk.ID", "Weight..g.")]



# Filter the relatedness matrix and euthanasia dataset based on common Chipmunk IDs
common_chipmunk_ids <- intersect(rownames(relatedness_matrix), euthanasia_subset$Chipmunk.ID)

# Subset the datasets to include only common individuals
filtered_euthanasia_data <- filtered_euthanasia_data[filtered_euthanasia_data$Chipmunk.ID %in% common_chipmunk_ids, ]
filtered_relatedness_matrix <- filtered_relatedness_matrix[common_chipmunk_ids, common_chipmunk_ids]
filtered_relatedness_matrix <- as.matrix(filtered_relatedness_matrix)

# Create a symmetric matrix with original diagonal values
symmetric_relatedness_matrix <- matrix(0, nrow = nrow(filtered_relatedness_matrix), ncol = ncol(filtered_relatedness_matrix))
upper_tri <- upper.tri(symmetric_relatedness_matrix)
lower_tri <- lower.tri(symmetric_relatedness_matrix)

symmetric_relatedness_matrix[upper_tri] <- filtered_relatedness_matrix[upper_tri]
symmetric_relatedness_matrix[lower_tri] <- t(filtered_relatedness_matrix)[lower_tri]

# Set diagonal values to the original diagonal values
diag(symmetric_relatedness_matrix) <- diag(filtered_relatedness_matrix)

# Reorder filtered_euthanasia_data
filtered_euthanasia_data <- filtered_euthanasia_data[match(common_chipmunk_ids, filtered_euthanasia_data$Chipmunk.ID), ]


# Print the filtered datasets
print(symmetric_relatedness_matrix)
print(filtered_euthanasia_data)

# Calculate the Euclidean distance matrix from weight differences
weight_diff_matrix <- as.matrix(dist(filtered_euthanasia_data$Weight..g., diag = TRUE, upper = TRUE))

# Run a Mantel test
mantel_test_result <- vegan::mantel(
  weight_diff_matrix,
  symmetric_relatedness_matrix,
  na.rm = TRUE
)



# Try this out with 1-relatedness
rel_minus1 <- 1 - symmetric_relatedness_matrix
mantel_test_result_Inverse_Relatedness <- vegan::mantel(
  weight_diff_matrix,
  rel_minus1,
  na.rm = TRUE
)


# Print the result of the Mantel test
print(mantel_test_result)
print(mantel_test_result_Inverse_Relatedness)


