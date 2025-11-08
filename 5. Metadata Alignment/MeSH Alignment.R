rm(list = ls(all.names = TRUE)) # clears global environ

# Load packages
library(tidyverse) # for data cleaning
library(dplyr) 

##########################################################################################################################
############################ Import sample metadata and remove unwanted diseases and tissues ############################# 
##########################################################################################################################

# Import data
load("Selected_Defined_Samples.RData")

# Remove samples with Null designation for tissue
Selected_Samples <- Selected_Samples[which(Selected_Samples$Tissue_Type != 'Null'), ]

# Isolate control samples from diseased samples
Control_Selected_Samples <- Selected_Samples[which(Selected_Samples$Disease_Status == 'Control'), ]
Selected_Samples <- Selected_Samples[which(Selected_Samples$Disease_Status == 'Disease'), ]

# Define desired diseases for isolation
Desired_Diseases <- c('Active Tuberculosis', 'Allergic Asthma', 'Allergic Rhinitis', 'Allergy', 'Aspergillosis', 
                      'Asthma', 'Bronchitis', 'Bronchopulmonary Dysplasia', 'Chronic Obstructive Pulmonary Disease', 
                      'Chronic Rhinosinusitis', 'Community Acquired Pneumonia', 'Cystic Fibrosis', 'Eosinophilic Asthma', 
                      'Idiopathic Pulmonary Fibrosis', 'Infant Acute Respiratory Distress Syndrome',  'Influenza', 
                      'Influenza A (H1N1)', 'Invasive Aspergillosis',  'Latent Tuberculosis', 'Lower Respiratory Tract Infections', 
                      'Lung Cancer', 'Lung Disease',  'Lung Squamous Cell Carcinoma', 'Lung Tumor',  'Mucormycosis', 'Non-Small Cell Lung Cancer', 
                      'Pneumonia', 'Pulmonary Fibrosis', 'Pulmonary Sarcoidosis', 'Rhinitis', 'Rhinosinusitis', 'Rhinovirus Infection', 
                      'Scleroderma-Associated Interstitial Lung Disease', 'Sepsis', 'Septic Shock', 'Small Cell Lung Cancer',
                      'Tracheal Stenosis', 'Tuberculosis', 'Tuberculosis Infection', 'Tuberculous Meningitis', 'Unspecified Septicemia')

# Select samples from desired disease list
Selected_Samples <- Selected_Samples[which(Selected_Samples$Disease_Name %in% Desired_Diseases), ]

# Remove food allergies
Selected_Samples <- Selected_Samples[which(Selected_Samples$Tissue_Type != 'Airway smooth muscle' & 
                                             Selected_Samples$Tissue_Type != 'Liver' & 
                                             Selected_Samples$Tissue_Type != 'LAD2 MCs cultured with Asthmatic AECs'), ]

# Select for control samples from tissues of selected samples
Control_Selected_Samples <- Control_Selected_Samples[which(tolower(Control_Selected_Samples$Tissue_Type) %in% tolower(Selected_Samples$Tissue_Type)), ]

##########################################################################################################################
########################### Import and format MeSH tree positions and corresponding MeSH terms ########################### 
##########################################################################################################################

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
                                                              'Tissue_Type',
                                                              MesH_Position_Number)]) %>% na.omit()
  
  # Isolate third tree fork
  Sample_Metadata[, 6] <- sub("^([^.]*.[^.]*.[^.]*).*", "\\1", Sample_Metadata[, 6])
  
  colnames(Sample_Metadata)[6] <- 'Disease_Designation'
  
  Sample_Metadata
  
})

# Combine every tree position for each sample to one spot
Sample_Metadata <- bind_rows(Sample_Metadata)

# Add tree terms for isolated tree fork
Sample_Metadata <- merge(Sample_Metadata, mtrees2024, by = 'Disease_Designation', relationship = "many-to-many")

# Rearrange columns
Sample_Metadata <- Sample_Metadata[, c('geo_accession', 
                                       'series_id', 
                                       'Tissue_Type', 
                                       'Disease_Name',
                                       'Used_MeSH_Term', 
                                       'Tree_Term')]

# Set to upper case to prevent capitalization from affecting grouping
Sample_Metadata <- sapply(Sample_Metadata, toupper)
Sample_Metadata <- Sample_Metadata %>% data.frame()

# Arrange samples by tissue type
Sample_Metadata <- Sample_Metadata %>% arrange(Tissue_Type)

########

# Set control tissue type as Tree Term and Disease_Designation
Control_Selected_Samples$Tree_Term <- Control_Selected_Samples$Disease_Status

# Rearrange columns to match diseased metadata
Control_Selected_Samples <- Control_Selected_Samples[, colnames(Sample_Metadata)]

# Set to upper case to prevent capitalization from affecting grouping
Control_Selected_Samples <- sapply(Control_Selected_Samples, toupper)
Control_Selected_Samples <- Control_Selected_Samples %>% data.frame()

# Arrange samples by tissue type
Control_Selected_Samples <- Control_Selected_Samples %>% arrange(Tissue_Type)

# Add tissue type information to tree term
Control_Selected_Samples$Tree_Term <- paste0(Control_Selected_Samples$Tree_Term, ' - ', Control_Selected_Samples$Tissue_Type)

# Remove samples with identical tree term designations
# As some tree terms will map to the same location based on tree split point selected
# This will prevent the same sample from appearing with the same labeling and improperly inflating tree term numbers
Sample_Metadata <- Sample_Metadata[!duplicated(paste0(Sample_Metadata$geo_accession, Sample_Metadata$Tree_Term)), ]

# Reset row names to be row number
rownames(Sample_Metadata) <- seq(1, nrow(Sample_Metadata), 1)

##########################################################################################################################
######################### Count samples by tree term and remove terms with less than 10 samples ########################## 
##########################################################################################################################

# Add control samples to dataset
Sample_Metadata <- rbind(Control_Selected_Samples, Sample_Metadata)

# Reset row names to be row number
rownames(Sample_Metadata) <- seq(1, nrow(Sample_Metadata), 1)

# Set tree term to a factor with set order to keep control samples as earlier categorical positions
Sample_Metadata$Tree_Term <- factor(Sample_Metadata$Tree_Term, levels = unique(Sample_Metadata$Tree_Term))

# Count number of samples in each polarization ID group
group_num <- Sample_Metadata %>% 
  group_by(Tree_Term) %>% 
  summarise(total_count = n(),
            .groups = 'drop')

# Filter out samples with less than 10 samples and remove unusable groups
group_num <- dplyr::filter(group_num, as.numeric(total_count) >= 10)

# Remove unwanted group
group_num <- group_num[which(group_num$Tree_Term != 'CENTRAL NERVOUS SYSTEM BACTERIAL INFECTIONS'), ]

# Remove filtered samples
Sample_Metadata <- Sample_Metadata[which(Sample_Metadata$Tree_Term %in% group_num$Tree_Term), ]

# Add numeric value representing each tree term. Python uses 0 indexing, so subtract 1.
Sample_Metadata$Tree_Term_Categorical <- as.numeric(factor(Sample_Metadata$Tree_Term)) - 1

save(Sample_Metadata, file = 'MeSH_aligned_metadata.RData')
