rm(list = ls(all.names = TRUE)) # clears global environ

# Load packages
library(tidyverse) # for data cleaning
library(dplyr) 
library(edgeR)

setwd("")

# Import data
load('Selected_Samples.RData')

setwd("")

# Disease names in DiSignDB do not necessarily align with MeSH terms.
# Terms must be manually selected based on key word matching and literature searches.
# write.csv(file = 'Disease_names.csv', unique(Selected_Samples$Disease_Name))

MeSH_Added <- read.csv('Disease_names - MeSH Added.csv')

# Change cells which did not properly import NA to NA
MeSH_Added[MeSH_Added == ""] <- NA

Selected_Samples <- merge(Selected_Samples, MeSH_Added, by = 'Disease_Name', relationship = 'many-to-many')
rownames(Selected_Samples) <- Selected_Samples$geo_accession

#############################################################

DiSignAtlas_Samples <- DiSignAtlas_Samples[, rownames(Selected_Samples)]

Selected_Samples$Disease_Status_Designation <- paste0(Selected_Samples$Disease_Status,
                                                      '_',
                                                      Selected_Samples$Used_MeSH_Term)

# Create DGEList object
d0 <- DGEList(DiSignAtlas_Samples,
              group = factor(Selected_Samples$Disease_Status_Designation))

# Calculate normalization factor
d0 <- calcNormFactors(d0)

# Use EdgeR function to automatically filter genes
# This gives logical values to keep or remove
keep.exprs <- filterByExpr(d0,
                           group = factor(Selected_Samples$Disease_Status_Designation))

# Remove filtered genes
d1 <- d0[keep.exprs, ]

# Calculate normalization factor
d1 <- calcNormFactors(d1)

# Calculate logCPM of filtered gene set
lcpm <- cpm(d1, log = TRUE)

save(lcpm, file = 'LCPM.RData')

save(Selected_Samples, file = 'Selected_Samples_Metadata.RData')