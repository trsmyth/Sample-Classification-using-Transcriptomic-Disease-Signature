rm(list = ls(all.names = TRUE)) # clears global environ.

# Load packages
library(tidyverse) # for data cleaning
library(dplyr) 
library(RSelenium)
library(rvest)

setwd("")

# Import sample designations from DiSignAtlas. 
# DiSignAtlas was accessed on 11/23/2024 including 
# samples linked to each MeSH term category, 
# sample metadata, and raw data files
Immune <- read.csv('Immune_System_Diseases.csv')
Infection <- read.csv('Infection.csv')
Nutrition <- read.csv('Nutritional_and_Metabolic.csv')
Respiratory <- read.csv('Respiratory_Tract_Diseases.csv')

# Compile all data to single list
All_sample_designations <- list(Immune, Infection, Nutrition, Respiratory)

# Add all data to single object
All_sample_designations <- do.call(rbind, All_sample_designations)

# Isolate Human RNAseq samples
All_sample_designations <- All_sample_designations[which(str_detect('Homo sapiens', All_sample_designations$Organism) &
                                                           str_detect('RNA-Seq', All_sample_designations$Data.Type)), ]

# Remove duplicated study IDs
All_sample_designations <- All_sample_designations[!duplicated(All_sample_designations$study.ID), ]

# Set row names to match study ID designation
rownames(All_sample_designations) <- All_sample_designations$study.ID

###################################################

# Import isolated sample metadata
Disease_information_Datasets <- read.csv('Disease_information_Datasets.csv')

# Isolate rows with study ID information which matches above isolated study IDs
Disease_information_Datasets <- Disease_information_Datasets[which(Disease_information_Datasets$dsaid %in% All_sample_designations$study.ID), ]

# Set row names to study ID
rownames(Disease_information_Datasets) <- Disease_information_Datasets$dsaid

# Rearrange rows to match between data frames
Disease_information_Datasets <- Disease_information_Datasets[rownames(All_sample_designations), ]

######################################################################################################
######################################################################################################
######################################################################################################

setwd("")

# Set up firefox as driver
rD <- rsDriver(browser = "firefox",
               chromever = NULL)

# Set up client
remDr <- rD$client

# Perform webscrape to isolate GEO IDs for control and diseased samples
Disease_Designation <- lapply(Disease_information_Datasets$dsaid, function(Single_ID){

  # Read html for all tr elements and isolate text
  html <- read_html(paste0("http://www.inbirg.com/disignatlas/detail/", Single_ID)) %>%
    html_nodes("tr") %>%
    html_text2()

  # Isolate GSM data
  Case_info <- html[which(str_detect(html, 'Control \\| Case'))]

  # Reformat control/disease data to arrange samples later
  Case_info <- Case_info %>%

    # Remove 'Control | Case\t'
    gsub(pattern = 'Control \\| Case\\\t',
         replacement = '') %>%

    # Remove '\t' in middle of string separating control from case box
    gsub(pattern = '\\\t.',
         replacement = ';G') %>%

    # Remove '\t' at end of string
    gsub(pattern = '\\\t',
         replacement = '') %>%

    gsub(pattern = 'Click to open',
         replacement = '') %>%

    gsub(pattern = '\\\n',
         replacement = '')

  # Split character sting into a list
  Case_info <- as.list(strsplit(Case_info, ';')) %>% unlist()

})

# Close the server
remDr$close()

names(Disease_Designation) <- Disease_information_Datasets$dsaid

save(Disease_Designation, file = 'Disease_Designation.RData')

######################################################################################################
######################################################################################################
######################################################################################################

load("Disease_Designation.RData")

# Add GSE data
All_sample_designations$GSE <- Disease_information_Datasets$accession

# Add number of control samples in dataset for each study id
All_sample_designations$Control_samples <- sub(pattern = '\\|.*', 
                                               replacement = '', 
                                               Disease_information_Datasets$control_case_sample_count)

# Add number of diseased samples in dataset for each study id
All_sample_designations$Disease_samples <- sub(pattern = '.*\\|', 
                                               replacement = '', 
                                               Disease_information_Datasets$control_case_sample_count)

# Remove unneeded data
rm(list = setdiff(ls(), c('All_sample_designations', 
                          'Disease_Designation')))

###################################################

# Create a list of the files from your target directory
file_list <- All_sample_designations$study.ID

# Initiate a blank data frame, each iteration of the loop will append the data from the given file to this variable
Samples <- list()

# For each file name, read in the corresponding csv file, skipping the first row as colnames are in row 2
for(i in 1:length(file_list)){
    
  # Import csv data
  Imported_data <- read.csv(paste0('', # Add directory here
                                   file_list[i], 
                                   '_profile.csv'), 
                            skip = 1)
  
  # Isolate the sample IDs - these are GSM designations
  Sample_ID <- colnames(Imported_data)[2:ncol(Imported_data)]
  
  # Define samples as control or diseased based on metadata description from DiSignAtlas
  Control_Disease <- c(rep('Control', 
                           as.numeric(All_sample_designations[which(All_sample_designations$study.ID == file_list[[i]]), 
                                                              'Control_samples'])), 
                       
                       rep('Disease', 
                           as.numeric(All_sample_designations[which(All_sample_designations$study.ID == file_list[[i]]), 
                                                              'Disease_samples'])))
  
  ifelse(length(Sample_ID) != length(Control_Disease),
         
         Samples[[i]] <- list(Sample_ID, Control_Disease),
         
         # Combine sample names, control/disease designations, sample disease name, and sample tissue type
         Samples[[i]] <- data.frame('Study_ID' = file_list[i],
                                    'Sample_ID' = Sample_ID[order(match(Sample_ID, Disease_Designation[[i]]))],
                                    'Disease_Status' = Control_Disease, 
                                    'Disease_Name' = rep(All_sample_designations[file_list[[i]], 'Disease.Name'], length(Sample_ID)), 
                                    'Tissue_Type' = rep(All_sample_designations[file_list[[i]], 'Tissue.Cell.Type'], length(Sample_ID))))
  
}

# Match sample names to file list
names(Samples) <- file_list

# Isolate data where number of samples imported does not match the number of samples listed in metadata
Unmatching <- Samples[sapply(Samples, class) != "data.frame"]

# Isolate data where number of samples DOES match number of samples listed in metadata
Samples <- Samples[sapply(Samples, class) == "data.frame"]

###################################################

# Combine data to single data frame
Samples_of_Interest <- do.call(rbind, Samples)

# Check if any samples are duplicated
Duplicated_samples <- Samples_of_Interest[duplicated(Samples_of_Interest$Sample_ID), ]

# Remove any samples which are duplicated - This keeps first instance of the duplicated sample
Samples_of_Interest <- Samples_of_Interest[!duplicated(Samples_of_Interest$Sample_ID), ]

save(Samples_of_Interest, file = 'Samples_of_Interest.RData')