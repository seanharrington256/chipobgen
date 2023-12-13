# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("LEA")
# install.packages(sequoia)
# install.packages(adegenet)

setwd("~/UWyo/kellyn_chipmunks")

# Load up packages
library(LEA)
library(sequoia)
library(adegenet)

# set up file paths
path_ugeno <- "chip_c93_branch1.u.geno" # unlinked geno file - I manually renamed it
path_ustr <- "chip_c93_branch1.ustr" # unlinked structure file
path_sequoia_file <- "chip_c93_branch1Sequoia.u.geno" # file to be made for Sequoia

# Read in the geno file
geno <- read.geno(path_ugeno)
nums_snps <- ncol(geno) # get the number of snps
num_ind <- nrow(geno) # get the number of individuals

# replace 9 with -9 across the matrix, then transpose it
geno_conv <- apply(geno, 1, function(x) as.numeric(gsub(9, -9 , x)))
tgeno_conv <- t(geno_conv)

  
# read in the structure file to get names of individuals
path_stru <- gsub(".ustr", ".stru", path_ustr)  # Use a regular expression substitution to generate the new file name
file.copy(path_ustr, path_stru) # make a copy of the file with the new name
# Now we can read in this file
ustr <- read.structure(path_stru, n.ind=num_ind, n.loc=nums_snps, onerowperind = FALSE, col.lab=1, col.pop=0, NA.char="-9", pop=NULL, ask=FALSE, quiet=FALSE)
ind_names <- rownames(ustr@tab) ## get the individual names in the order that they show up in the various files - this is important farther down for getting coordinates into the right order for plotting

# add individual names as row names
rownames(tgeno_conv) <- ind_names
# write that as a table
write.table(tgeno_conv, path_sequoia_file, col.names = FALSE, quote = FALSE)


# Read it into Sequoia:
GenoM <- as.matrix(read.table(path_sequoia_file, row.names=1, header=FALSE))



write.table(ind_names, "retained_individuals.csv", col.names = FALSE, row.names = FALSE)



