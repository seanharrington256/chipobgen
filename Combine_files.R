
# setwd("C:/Users/UW-User/Desktop/")
setwd("~/UWyo/kellyn_chipmunks/Run_sequoia") # Sean's working dir
library(sequoia)
library(kinship2)
library(dplyr)

#setting up file paths
path_sex_birth_file <- "Sequoia_Sex_Birth_Data_1504_fixed.csv"
path_seq_out <- "pedigree.csv"
loci <- "Sample_cov.csv"
euth_log <- "Euthanasia_Log.csv"

# read in the files
ped_info <- read.csv(path_seq_out, row.names = 1)
sex_birth <- read.csv(path_sex_birth_file, row.names = 1)
coverage <- read.csv(loci)
euth <- read.csv(euth_log)
weights <- dplyr::select(euth, c(Chipmunk.ID, Weight..g., Enclosure))
colnames(weights) <- c("ID", "weight_g", "encolsure")
weights$weight_g <- as.numeric(gsub("N/A", NA, weights$weight_g))


# merge them
ped_sexbirth <- merge(ped_info, sex_birth, by.x = "id", by.y = "ID", all = TRUE)

# individuals in sex_birth not in ped_info and the reverse
ped_sexbirth$id[!ped_sexbirth$id %in% sex_birth$ID]
sex_birth$ID[!sex_birth$ID %in% ped_sexbirth$id]


# merge with coverage
ped_sexbirth_loci <- merge(ped_sexbirth, coverage, by.x = "id", by.y = "ID", all = TRUE)

#find any non-overlapping IDs
ped_sexbirth$id[!ped_sexbirth$id %in% coverage$ID]
coverage$ID[!coverage$ID  %in% ped_sexbirth$id]

# merge with weights

# make some replacements
weights$ID <- gsub("untagged female", "untagged_female", weights$ID)

in_ped_not_weight <- ped_sexbirth_loci$id[!ped_sexbirth_loci$id %in% weights$ID]
in_weight_not_ped <- weights$ID[!weights$ID %in% ped_sexbirth_loci$id]

all_data <- merge(ped_sexbirth_loci, weights, by.x = "id", by.y = "ID", all = TRUE)


in_ped_not_weight <- ped_sexbirth_loci$id[!ped_sexbirth_loci$id %in% weights$ID]
in_weight_not_ped <- weights$ID[!weights$ID %in% ped_sexbirth_loci$id]


write.csv(all_data, file = "chipmunk_loci_etc.csv")
