rm(list = ls(all.names = TRUE)) # clears global environ

# Load packages
library(tidyverse) # for data cleaning
library(dplyr) 
library(edgeR)

setwd("")

# Import data
load("Selected_Defined_Samples.RData")

# Change Disease_Status_Designation to list control tissue type or disease name
Selected_Samples <- Selected_Samples %>% mutate(Disease_Status_Designation = 
                                                  case_when(str_detect(Disease_Status, 'Control') ~ paste0('Control ', Tissue_Type), 
                                                
                                                TRUE ~ as.character(.$Used_MeSH_Term)))

##########################################################################################################################
########################### Import and format MeSH tree positions and corresponding MeSH terms ########################### 
##########################################################################################################################

setwd("")

# Convert raw xml data downloaded from MeSH to an r data frame.
# This can be drastically sped up if xml structural information added to function.
# supp2024 <- xmlToDataFrame("supp2024.xml")
# desc2024 <- xmlToDataFrame("desc2024.xml")

# Import MeSH tree terms and positions
mtrees2024 <- read.csv('mtrees2024.csv', header = FALSE)

# Isolate tree terms and positions
Tree_term <- sub(pattern = ';.*', '', mtrees2024$V1)
Tree_position <- sub(pattern = '.*;', '', mtrees2024$V1)

# Combine tree terms and positions into single data frame
mtrees2024 <- data.frame(Tree_Term = Tree_term,
                         Disease_Designation = Tree_position)

##########################################################################################################################
############################## Format metadata by separating diseased from control samples ############################### 
########################### For each assigned tree position, isolate the indicated fork number ###########################
################################## and add the corresponding tree term to the metadata. ################################## 
################################ For control data, assign the tissue type as the tree term ############################### 
##########################################################################################################################

# Isolate single set of tree positions and cut tree positions to indicated fork
Sample_Metadata <- lapply(colnames(Selected_Samples)[which(str_detect(colnames(Selected_Samples), 'Tree_Position'))], function(MesH_Position_Number){
  
  Diseased_Selected_Samples <- Selected_Samples[which(Selected_Samples$Disease_Status == 'Disease'), ]

  # Isolate a single set of tree positions and remove NA results
  Sample_Metadata <- data.frame(Diseased_Selected_Samples[, c('Disease_Name', 
                                                              'Used_MeSH_Term', 
                                                              'geo_accession', 
                                                              'series_id',
                                                              'Disease_Status_Designation',
                                                              'Tissue_Type',
                                                              MesH_Position_Number)]) %>% na.omit()
  
  # Isolate fourth tree fork
  # Sample_Metadata[, 7] <- sub("^([^.]*.[^.]*.[^.]*.[^.]*).*", "\\1", Sample_Metadata[, 7])
  
  # Isolate third tree fork
  Sample_Metadata[, 7] <- sub("^([^.]*.[^.]*.[^.]*).*", "\\1", Sample_Metadata[, 7])
  
  # Isolate second tree fork
  # Sample_Metadata[, 7] <- sub("^([^.]*.[^.]*).*", "\\1", Sample_Metadata[, 7])
  
  colnames(Sample_Metadata)[7] <- 'Disease_Designation'
  
  Sample_Metadata
  
})

# Combine every tree position for each sample to one spot
Sample_Metadata <- do.call(rbind, Sample_Metadata)

# Add tree terms for isolated tree fork
Sample_Metadata <- merge(Sample_Metadata, mtrees2024, by = 'Disease_Designation', relationship = "many-to-many")

# Rearrange columns
Sample_Metadata <- Sample_Metadata[, c('geo_accession', 
                                       'series_id', 
                                       'Tissue_Type', 
                                       'Disease_Name',
                                       'Used_MeSH_Term', 
                                       'Tree_Term')]

# Set tree term designation to include tissue source
# This will mean the future model attempts to ID asthma - blood and asthma - airway epithelial cell separately
# This is an alternative attempt to batch effect removal (treating tissue as batch) which did not appear to work
Sample_Metadata$Tree_Term <- paste0(Sample_Metadata$Tree_Term, ' - ', Sample_Metadata$Tissue_Type)

# Isolate control samples metadata
Control_Selected_Samples <- Selected_Samples[which(Selected_Samples$Disease_Status == 'Control'), ]

# Set control tissue type as Tree Term and Disease_Designation
Control_Selected_Samples$Tree_Term <- Control_Selected_Samples$Disease_Status_Designation

# Rearrange columns to match diseased metadata
Control_Selected_Samples <- Control_Selected_Samples[, colnames(Sample_Metadata)]

# Combine metadata to one location
Sample_Metadata <- rbind(Sample_Metadata, Control_Selected_Samples)

# Remove 'Null' tissue designation
# This is caused by missing metadata
# Manual curation may be able to find tissue information
Sample_Metadata <- Sample_Metadata[which(Sample_Metadata$Tissue_Type != 'Null'), ]

# Reset row names to be row number
rownames(Sample_Metadata) <- seq(1, nrow(Sample_Metadata), 1)

# Count number of samples in each polarization ID group
group_num <- Sample_Metadata %>% 
  group_by(Tree_Term) %>% 
  summarise(total_count = n(),
            .groups = 'drop')

# Filter out samples with less than 10 samples and remove unusable groups
group_num <- dplyr::filter(group_num, as.numeric(total_count) >= 10)

##########################################################################################################################

##########################################################################################################################

load("")

# Transcribe the count data and convert to df
DiSignAtlas_Samples <- t(DiSignAtlas_Samples) %>%
  data.frame() %>%
  rownames_to_column()

colnames(DiSignAtlas_Samples)[1] <- 'geo_accession'

# Merge metadata with count data
DiSignAtlas_Samples <- merge(Sample_Metadata, DiSignAtlas_Samples, by = 'geo_accession', relationship = 'many_to_many')

# Create DGEList object
d0 <- DGEList(t(DiSignAtlas_Samples[, 7:ncol(DiSignAtlas_Samples)]),
              group = factor(DiSignAtlas_Samples$Tree_Term))

# Calculate normalization factor
d0 <- calcNormFactors(d0)

# Use EdgeR function to automatically filter genes
# This gives logical values to keep or remove
keep.exprs <- filterByExpr(d0,
                           group = factor(DiSignAtlas_Samples$Tree_Term))

# Remove filtered genes
d1 <- d0[keep.exprs, ]

# Calculate normalization factor
d1 <- calcNormFactors(d1)

# Calculate logCPM of filtered gene set
lcpm <- cpm(d1, log = TRUE)

# Create model matrix to remove batch effects
Design <- model.matrix(~ Tree_Term,
                       data = DiSignAtlas_Samples)

# Remove batch effects caused by experimental series.
# Requires samples as columns
Batched_lcpm_data <- limma::removeBatchEffect(t(lcpm),
                                              batch = DiSignAtlas_Samples$series_id,
                                              design = Design)

# Transcribe data and convert to df so adding metadata is faster
Batched_lcpm_data <- t(Batched_lcpm_data) %>% data.frame()

# gc to clear memory
gc()

# Add geo accession data
Batched_lcpm_data$geo_accession <- Sample_Metadata$geo_accession

# Move geo accession to first column
Batched_lcpm_data <- Batched_lcpm_data %>% select(geo_accession, everything())

# Reset row names to be row number
rownames(Batched_lcpm_data) <- seq(1, nrow(Batched_lcpm_data), 1)

setwd("")

save(Batched_lcpm_data, file = 'Batched_lcpm_data.RData')

# load('Batched_lcpm_data.RData')

##########################################################################################################################

##########################################################################################################################

# Add numeric value representing each tree term
# Python uses 0 indexing, so subtract 1
Sample_Metadata$Tree_Term_Categorical <- as.numeric(factor(Sample_Metadata$Tree_Term)) - 1

# Count number of samples in each polarization ID group
group_num <- Sample_Metadata %>% 
  group_by(Used_MeSH_Term) %>% 
  summarise(total_count=n(),
            .groups = 'drop')

# Filter out samples with less than 5 samples and remove unusable groups
group_num <- dplyr::filter(group_num, as.numeric(total_count) >= 5)

# Set up data splitting for each tree term
Splitting <- lapply(group_num$Used_MeSH_Term, function(Single_Term){
  
  # Isolate single tree term
  Sample_Metadata_Subset <- Sample_Metadata[which(Sample_Metadata$Sample_Metadata_TOI == Single_Term), ]
  
  # Split data into 70/30 split
  split <- caTools::sample.split(seq_len(nrow(Sample_Metadata_Subset)), SplitRatio = 0.7)
  
  # Create train and test data sets
  Training_Data <- subset(Sample_Metadata_Subset, split == 'TRUE')
  Testing_Data <- subset(Sample_Metadata_Subset, split == 'FALSE')
  
  list(Training_Data, 
       Testing_Data)
  
})

# Bind training and testing data together
Sample_Metadata <- do.call(Map, c(f = rbind, Splitting))

##########################################################################################################################

##########################################################################################################################

# Create training data by merging selected training data metadata with matching count data
Training_Data <- merge(Sample_Metadata[1], Batched_lcpm_data, by = 0)

# Create testing data by merging selected testing data metadata with matching count data
Testing_Data <- merge(Sample_Metadata[2], Batched_lcpm_data, by = 0)

#######

save(Training_Data, 
     Testing_Data, 
     file = 'LCPM_Data.RData')
