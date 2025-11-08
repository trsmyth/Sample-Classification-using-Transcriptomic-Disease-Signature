rm(list = ls(all.names = TRUE)) # clears global environ

# Load packages
library(tidyverse) # for data cleaning
library(dplyr) 
library(edgeR)

load("MeSH_aligned_metadata.RData")
load("DiSignAtlas_Data_ENSEMBL.RData")

DiSignAtlas_Samples <- t(DiSignAtlas_Samples) %>% data.frame() %>% rownames_to_column()
colnames(DiSignAtlas_Samples)[1] <- 'geo_accession'

# Merge metadata with count data
DiSignAtlas_Samples <- merge(Sample_Metadata, DiSignAtlas_Samples, by = 'geo_accession', relationship = 'many_to_many')

# Create a metadata data frame
Sample_Designation <- DiSignAtlas_Samples[, 1:7]
DiSignAtlas_Samples <- t(DiSignAtlas_Samples[, 8:ncol(DiSignAtlas_Samples)])

# Create DGEList object
d0 <- DGEList(DiSignAtlas_Samples,
              group = factor(Sample_Designation$Tree_Term))

# Calculate normalization factor
d0 <- calcNormFactors(d0)

# Use EdgeR function to automatically filter genes
# This gives logical values to keep or remove
keep.exprs <- filterByExpr(d0,
                           group = factor(Sample_Designation$Tree_Term))

# Remove filtered genes
d1 <- d0[keep.exprs, ]

# Calculate normalization factor
d1 <- calcNormFactors(d1)

# Calculate logCPM of filtered gene set
lcpm <- cpm(d1, log = TRUE)

lcpm <- lcpm %>% data.frame()

colnames(lcpm) <- seq(1, ncol(lcpm), 1)

########

# Isolate disease sample count and metadata
Disease_Sample_Designation <- Sample_Designation[which(!grepl("CONTROL", Sample_Designation$Tree_Term)), ]
Diseased <- lcpm[, rownames(Disease_Sample_Designation)]

# Create model matrix to remove batch effects
Design <- model.matrix(~ Tree_Term,
                       data = Disease_Sample_Designation)

# Remove batch effects caused by experimental series and tissue type to get universal 
# tree term expression profile. Requires samples as columns.
Diseased <- limma::removeBatchEffect(Diseased,
                                     batch = Disease_Sample_Designation$series_id,
                                     batch2 = Disease_Sample_Designation$Tissue_Type,
                                     design = Design)

########

# Isolate control sample count and metadata
Control_Sample_Designation <- Sample_Designation[which(grepl("CONTROL", Sample_Designation$Tree_Term)), ]
Control <- lcpm[, rownames(Control_Sample_Designation)]

# Create model matrix to remove batch effects
Design <- model.matrix(~ Tree_Term,
                       data = Control_Sample_Designation)

# Remove batch effects caused by experimental series. 
# Since control sample tree term requires tissue type, only correct for sample ID.
# Requires samples as columns.
Control <- limma::removeBatchEffect(Control,
                                    batch = Control_Sample_Designation$series_id,
                                    design = Design)

# Combine Diseased and Control batch corrected data and metadata
lcpm <- merge(Diseased, Control, by = 0)
lcpm <- data.frame(lcpm[-1], row.names = lcpm$Row.names)

colnames(lcpm) <- gsub(pattern = 'X', replacement = '', colnames(lcpm))

Sample_Designation <- Sample_Designation[colnames(lcpm), ]

##########################################################################################################################
#################################### Remove genes which demonstrate high within group, ################################### 
############################### low between group, and low between group:within group ratio ############################## 
##########################################################################################################################

library(collapse)

lcpm <- t(lcpm) %>% data.frame()

# Add tree term to batched data
lcpm <- data.frame(factor(Sample_Designation$Tree_Term), lcpm) %>% 
  setNames(c('Tree_Term', colnames(lcpm)))

# Calculate variance after removing batch effects to estimate effect of disease
Gene_var <- lcpm %>% 
  fgroup_by(Tree_Term) %>% 
  fsummarise(across(colnames(lcpm)[2:ncol(lcpm)], fvar))

# Calculate mean of variances previously calculated for each group
column_means <- colMeans(Gene_var[2:ncol(Gene_var)])

# Calculate the overall variance of lcpm data for each gene
column_var <- sapply(lcpm[2:ncol(lcpm)], var)

# Combine results into single data frame
Variance_data <- data.frame(Mean_var = column_means, 
                            Overall_var = column_var)

# Calculate the ratio of difference to mean variance
Variance_data$Difference_ratio <- Variance_data$Overall_var / Variance_data$Mean_var

############

# Remove genes with mean variance less than 5
Variance_data <- Variance_data[which(Variance_data$Mean_var < 5), ]

# Remove genes with difference ratio above 5 (i.e. 5x more between than within variance)
Variance_data <- Variance_data[which(Variance_data$Difference_ratio > 5), ]

# Retain genes identified based on within group variance, between group variance, and ratio between variances
lcpm <- lcpm[, rownames(Variance_data)]

Genes = colnames(lcpm)

##########################################################################################################################
################################### Perform dendrogram clustering for each tree term or ################################## 
################################ tree term and tissue to remove poorly correlating samples ###############################
##########################################################################################################################

dend <- lapply(unique(Sample_Designation$Tree_Term), function(Single_Term){
  
  tmp_metadata <- Sample_Designation[which(Sample_Designation$Tree_Term == Single_Term), ]
  tmp <- lcpm[rownames(tmp_metadata), ]
  names <- rownames(tmp)
  
  # Calculate sample correlation and create dendrograms
  tmp_dend <- cor(t(tmp))
  tmp_dend <- hclust(as.dist(1 - tmp_dend))

  # Create clusters
  clu = cutree(tmp_dend, 
               h = 0.5) # Height for cut (1 - correlation coefficient)
  
  # Add cluster number information to data
  tmp <- cbind(clu, tmp) %>% data.frame()
  colnames(tmp)[1] <- 'cluster'
  
  # Determine the largest cluster and isolate samples in that cluster
  largest_cluster = names(rev(sort(table(tmp$cluster))))[1]
  tmp[which(tmp$cluster == largest_cluster), ]
  
})

lcpm <- bind_rows(dend)

Sample_Designation <- Sample_Designation[rownames(lcpm), ]

lcpm <- lcpm[-1]

##########################################################################################################################

# Set up data splitting for each tree term
Splitting <- lapply(unique(Sample_Designation$Tree_Term), function(Single_Term){
  
  # Isolate single tree term
  Sample_Metadata_Subset <- Sample_Designation[which(Sample_Designation$Tree_Term == Single_Term), ]
  
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

# Create training data by merging selected training data metadata with matching count data
Training_Data <- merge(Sample_Metadata[1], lcpm, by = 0)

# Create testing data by merging selected testing data metadata with matching count data
Testing_Data <- merge(Sample_Metadata[2], lcpm, by = 0)

#######

save(Training_Data,
     Testing_Data, 
     file = 'Pulmonary_LCPM_Data.RData')