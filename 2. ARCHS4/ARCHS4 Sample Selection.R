rm(list = ls(all.names = TRUE)) # clears global environ

library(rhdf5)
library(dplyr)
library(stringr)
library(tidyverse)
library(tibble)

# Load data corresponding to sample of interest
load("")

colnames(Samples_of_Interest)[2] <- 'geo_accession'
rownames(Samples_of_Interest) <- Samples_of_Interest$geo_accession

destination_file = "" # Originally accessed on 11/23/2024
extracted_expression_file = "DiSignAtlas_Sample_Selected_expression_matrix.tsv"

# Check h5 file structure to ID extractable information
# h5dump(destination_file, load = FALSE)

# Retrieve information from compressed data
Samples = as.data.frame(h5read(destination_file, '/meta/samples'))

# Set row names to GSM ID
row.names(Samples) = Samples$geo_accession

# Create not in function by negating in function
`%ni%` <- negate(`%in%`)

# Check for samples which are not in both data sets
Unmatch <- Samples_of_Interest[which(Samples_of_Interest$geo_accession %ni% Samples$geo_accession), ]

# Isolate selected samples and key metadata columns
Selected_Samples <- Samples[Samples_of_Interest$geo_accession, 
                            c('characteristics_ch1', 
                              'geo_accession', 
                              'series_id', 
                              'source_name_ch1', 
                              'title')] %>% na.omit()

Selected_Samples <- merge(Selected_Samples, 
                          Samples_of_Interest[rownames(Selected_Samples), ], 
                          by = 'geo_accession')

# Remove unneeded data
rm(list = setdiff(ls(), c('Selected_Samples', 
                          'destination_file', 
                          'extracted_expression_file')))

######################################################

# Selected samples to be extracted
samp = as.character(Selected_Samples$geo_accession)

# Retrieve information from compressed data
samples = h5read(destination_file, "meta/samples/geo_accession")
genes = h5read(destination_file, "meta/genes/ensembl_gene")

# Identify columns to be extracted
sample_locations = which(samples %in% samp)

# Subset sample locations to feed into a loop for extraction
sample_location_subsets <- list(sample_locations[1:2500], 
                                sample_locations[2501:5000], 
                                sample_locations[5001:7500], 
                                sample_locations[7501:10000], 
                                sample_locations[10001:12500], 
                                sample_locations[12501:length(sample_locations)])

setwd("")

a <- Sys.time()

# Loop through sample location subsets to extract only a subset of samples at a time, writing each subset to its own file.
# Without loop, a 64 GB system was maxed out and crashed. In current form, peaks are around 35 GB.
# This is still an extremely slow process.
for(i in 1:length(sample_location_subsets)){
  
  # extract gene expression from compressed data
  expression = t(h5read(destination_file, 
                        "data/expression", 
                        index = list(sample_location_subsets[[i]], 
                                     1:length(genes))))
  
  rownames(expression) = genes
  colnames(expression) = samples[sample_location_subsets[[i]]]
  
  # Print file
  write.table(expression, 
              file = paste0(i, '_', extracted_expression_file), 
              sep = "\t", 
              quote = FALSE, 
              col.names = NA)
  
  print(paste0("Expression file was created at ", 
               getwd(), 
               "/", 
               paste0(i, '_', extracted_expression_file)))
  
}

H5close()

b <- Sys.time()
b-a

# Remove unneeded data
rm(list = setdiff(ls(), c('Selected_Samples')))

######################################################

# Import count data subsets
DiSignAtlas_Samples <- list(Subset_1 = read_tsv('1_DiSignAtlas_Sample_Selected_expression_matrix.tsv'),
                            Subset_2 = read_tsv('2_DiSignAtlas_Sample_Selected_expression_matrix.tsv'),
                            Subset_3 = read_tsv('3_DiSignAtlas_Sample_Selected_expression_matrix.tsv'),
                            Subset_4 = read_tsv('4_DiSignAtlas_Sample_Selected_expression_matrix.tsv'),
                            Subset_5 = read_tsv('5_DiSignAtlas_Sample_Selected_expression_matrix.tsv'), 
                            Subset_6 = read_tsv('6_DiSignAtlas_Sample_Selected_expression_matrix.tsv'))

# Combine count data subsets into single location
# Count data should be in gene(row) x sample(column) layout for future analysis
DiSignAtlas_Samples <- Reduce(function(x, y) merge(x, y, all = TRUE), DiSignAtlas_Samples)

# Rename columns to match sample names
colnames(DiSignAtlas_Samples) <- c('GeneID', 
                                   sub(pattern = '.*\\.', 
                                       replacement = '', 
                                       colnames(DiSignAtlas_Samples)[2:ncol(DiSignAtlas_Samples)]))

# Set rownames to GeneID
rownames(DiSignAtlas_Samples) <- DiSignAtlas_Samples$GeneID

# Remove GeneID column
DiSignAtlas_Samples <- DiSignAtlas_Samples[, -1]

# Set rownames to GSM ID
rownames(Selected_Samples) <- Selected_Samples$geo_accession

save(DiSignAtlas_Samples, 
     file = 'DiSignAtlas_Data_ENSEMBL.RData')

save(Selected_Samples,
     file = 'Selected_Samples.RData')
