rm(list = ls(all.names = TRUE)) # clears global environ

# Load packages
library(tidyverse) # for data cleaning
library(dplyr) 
library(edgeR)

# Import data
load('Selected_Samples.RData')

# Disease names in DiSignDB do not necessarily align with MeSH terms.
# Terms must be manually selected based on key word matching and literature searches.
write.csv(file = 'Disease_names.csv', unique(Selected_Samples$Disease_Name))

MeSH_Added <- read.csv('Disease_names - MeSH Added.csv')

# Change cells which did not properly import NA to NA
MeSH_Added[MeSH_Added == ""] <- NA

Selected_Samples <- merge(Selected_Samples, MeSH_Added, by = 'Disease_Name', relationship = 'many-to-many')
rownames(Selected_Samples) <- Selected_Samples$geo_accession

save(Selected_Samples, file = 'Selected_Samples_Metadata.RData')