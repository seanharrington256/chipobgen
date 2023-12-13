
# setwd("C:/Users/UW-User/Desktop/")
setwd("~/UWyo/kellyn_chipmunks/Run_sequoia") # Sean's working dir
library(sequoia)
library(kinship2)
library(dplyr)

#setting up file paths
path_sequoia_file <- "chip_c93_branch1Sequoia.u.geno" # file to be made for Sequoia
path_sex_birth_file <- "Sequoia_Sex_Birth_Data_1504_fixed.csv"

#read in genetic data
GenoM <- as.matrix(read.table(path_sequoia_file, row.names=1, header=FALSE))
Sex_Birth <- read.csv(path_sex_birth_file)[,-1]

#  duplicate check & parentage assignment (takes few minutes)
ParOUT <- sequoia(GenoM = GenoM,  LifeHistData = Sex_Birth, Err=0.005,
                  Module="par", quiet = FALSE, Plot = TRUE)

ParOUT$DupGenotype  # inspect duplicates (intentional or accidental)

# polish dataset: 
# remove duplicate indiv. 985113005719233
Geno2 <- GenoM[!rownames(GenoM) %in% "985113005719233",]
Sex_Birth2 <- Sex_Birth[!Sex_Birth$ID %in% "985113005719233",]

# & drop low call rate samples
# & drop SNPs low MAF
stats <- SnpStats(GenoM, ParOUT$PedigreePar)
MAF <- ifelse(stats[,"AF"] <= 0.5, stats[,"AF"], 1-stats[,"AF"])
Geno2 <- Geno2[, -which(MAF < 0.05)]
#
# Indiv.Mis <- apply(Geno2, 1, function(x) sum(x == -9)) / ncol(Geno2)
# Geno2 <- Geno2[Indiv.Mis < 0.2, ]

### Second duplicate check & parentage assignment (takes few minutes)
ParOUT2 <- sequoia(GenoM = Geno2,  LifeHistData = Sex_Birth2, Err=0.005,
                  Module="par", quiet = FALSE, Plot = TRUE)

ParOUT2$DupGenotype  # inspect duplicates (intentional or accidental)



#run sequoia!!!!
SeqOUT <-sequoia(GenoM = Geno2, LifeHistData = Sex_Birth2, Err=0.005)


#save output
save(SeqOUT, file="Sequoia_output_10_25.RData")

SummarySeq(SeqOUT)


# ## Try removing individuals with more than 23% missing data
# Indiv.Mis <- apply(Geno2, 1, function(x) sum(x == -9)) / ncol(Geno2)
# Geno3 <- Geno2[Indiv.Mis < 0.23, ]
# 
# ParOUT3 <- sequoia(GenoM = Geno3,  LifeHistData = Sex_Birth2, Err=0.005,
#                    Module="par", quiet = FALSE, Plot = TRUE)
# 
# # run sequoia
# SeqOUT2 <-sequoia(GenoM = Geno3, LifeHistData = Sex_Birth2, Err=0.005)
# SummarySeq(SeqOUT2)
#   ## Doesn't really seem better to me??


# Try removing SNPs with high missing data
Geno4 <- GenoM[!rownames(GenoM) %in% "985113005719233",]

# & drop low call rate samples
# & drop SNPs low MAF
stats <- SnpStats(Geno4, ParOUT$PedigreePar)
MAF <- ifelse(stats[,"AF"] <= 0.5, stats[,"AF"], 1-stats[,"AF"])
Geno4 <- Geno4[, -which(MAF < 0.05 | stats$Mis > 0.2)]

ParOUT2 <- sequoia(GenoM = Geno4,  LifeHistData = Sex_Birth2, Err=0.005,
                   Module="par", quiet = FALSE, Plot = TRUE)

Indiv.Mis4 <- apply(Geno4, 1, function(x) sum(x == -9)) / ncol(Geno4)

stats2 <- SnpStats(Geno4, ParOUT2$PedigreePar)



SeqOUT3 <-sequoia(GenoM = Geno4, LifeHistData = Sex_Birth2, Err=0.005)
SummarySeq(SeqOUT3)

save(SeqOUT3, file="Sequoia_output_filtered.RData")


SeqOUT3$Pedigree


# make an object for kinship2 to plot out a pedigree
sex_birth_red <- Sex_Birth[match(SeqOUT3$Pedigree$id, Sex_Birth$ID), ] # match up sexbirth with the pedigree
plot_data1 <- merge(sex_birth_red, SeqOUT3$Pedigree, by.x = "ID", by.y = "id")[, c("ID", "Sex", "BirthYear", "dam", "sire")] # merge and keep only relevant columns
# add a column to dummy individuals for birth year (think this is unnecessary, but might be good to look at)
dummys <- SeqOUT3$DummyIDs
dummys$BirthYear <- NA
# and only keep relevant columns
dummys <- dummys[, c("id", "Sex", "BirthYear", "dam", "sire")]
colnames(dummys) <- c("ID", "Sex", "BirthYear", "dam", "sire")
# combine
plot_data <- rbind(plot_data1, dummys)

  
plot_data$Sex[plot_data$Sex == 2] <- "male"
plot_data$Sex[plot_data$Sex == 1] <- "female"



# process/reduce this so that if individuals have either a mother or father, but not both,
#   that both get assigned NA
plot_data_red <- plot_data %>%
  rowwise() %>%
  mutate(across(c(dam, sire), ~ ifelse(any(is.na(c(dam, sire))), NA, .)))

# untagged_female to untagF
plot_data_red$ID <- gsub("untagged_female", "untagF", plot_data_red$ID)
plot_data_red$dam <- gsub("untagged_female", "untagF", plot_data_red$dam)
plot_data_red$sire <- gsub("untagged_female", "untagF", plot_data_red$sire)
# Get last 4 digits of long IDs starting with 98
to_pull4 <- plot_data_red$ID[grep("^98", plot_data_red$ID)]
plot_data_red$ID[grep("^98", plot_data_red$ID)] <- substr(to_pull4, nchar(to_pull4) - 3, nchar(to_pull4))

# Function to apply the transformation
transform_ids <- function(id_column) {
  to_pull4 <- id_column[grep("^98", id_column)]
  id_column[grep("^98", id_column)] <- substr(to_pull4, nchar(to_pull4) - 3, nchar(to_pull4))
  return(id_column)
}

# Apply the transformation to specific columns
plot_data_red <- plot_data_red %>%
  mutate(across(c("ID", "dam", "sire"), transform_ids))


pedAll <- pedigree(id=plot_data_red$ID, 
                   dadid=plot_data_red$sire, momid=plot_data_red$dam,
                   sex=plot_data_red$Sex)

pdf(file = "chip_pedigree.pdf", width = 50, height = 5)
plot(pedAll)
dev.off()

write.csv(SeqOUT3$Pedigree, file = "pedigree.csv")

