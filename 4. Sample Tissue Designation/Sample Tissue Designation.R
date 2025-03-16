rm(list = ls(all.names = TRUE)) # clears global environ

# Load packages
library(tidyverse) # for data cleaning
library(dplyr) 

setwd("")

# Import data
load("Selected_Samples_Metadata.RData")

Selected_Samples$Filtered_Tissue_Type <- NA

Selected_Samples <- Selected_Samples %>% select(c(source_name_ch1, Tissue_Type, characteristics_ch1, Used_MeSH_Term, Filtered_Tissue_Type), everything())

#########################################

Selected_Samples <- Selected_Samples %>% mutate(Tissue_Type = 
                                                  case_when(str_detect(toupper(Tissue_Type), 'WHOLE BLOOD') |
                                                              str_detect(toupper(Tissue_Type), 'PERIPHERAL BLOOD') | 
                                                              str_detect(toupper(Tissue_Type), 'VENOUS BLOOD') ~ 'Blood',
                                                    
                                                    TRUE ~ as.character(.$Tissue_Type)))

#########################################

Selected_Samples <- Selected_Samples %>% mutate(Filtered_Tissue_Type =
                                                case_when(str_detect(toupper(source_name_ch1), 'FIBROBLAST') ~ 'Fibroblast',
                                                          
                                                          str_detect(toupper(source_name_ch1), ' T CELL') |
                                                            str_detect(toupper(source_name_ch1), ' T-CELL') |
                                                            str_detect(toupper(source_name_ch1), 'CD8 T') |
                                                            str_detect(toupper(source_name_ch1), 'TCELL: CD4') |
                                                            str_detect(toupper(source_name_ch1), 'CELL TYPE: CD8+ T') |
                                                            str_detect(toupper(source_name_ch1), 'TH CELL') |
                                                            str_detect(toupper(source_name_ch1), 'TSCM CELLS') |
                                                            str_detect(toupper(source_name_ch1), 'TH2 CELLS') |
                                                            str_detect(toupper(source_name_ch1), 'JURKAT') |
                                                            str_detect(toupper(source_name_ch1), 'EFFECTOR MEMORY TH1 CELLS') |
                                                            str_detect(toupper(source_name_ch1), 'EFFECTOR MEMORY TREG CELLS') |
                                                            str_detect(toupper(source_name_ch1), 'CULTURED T LYMPHOCYTES') |
                                                            str_detect(characteristics_ch1, 'TISSUE: BLOOD, CELL TYPE: CD3+CD4+CD45RA-CCR7+CD25-CCR4+ T CELLS') ~ 'T cells',
                                                          
                                                          str_detect(toupper(source_name_ch1), 'B CELL') | 
                                                            str_detect(toupper(source_name_ch1), 'CD19') | 
                                                            str_detect(toupper(source_name_ch1), 'B-CELL') ~ 'B cells',
                                                          
                                                          str_detect(toupper(source_name_ch1), 'WHOLE BLOOD LEUKOCYTE') |
                                                            str_detect(toupper(source_name_ch1), 'WHITE BLOOD CELL') ~ 'Leukocyte',
                                                          
                                                          str_detect(toupper(source_name_ch1), 'VENOUS BLOOD') ~ 'Eosinophil',
                                                          
                                                          str_detect(toupper(source_name_ch1), 'VASTUS LATERALIS') |
                                                            str_detect(toupper(source_name_ch1), 'MYOBUNDLE') ~ 'Muscle',
                                                          
                                                          str_detect(toupper(source_name_ch1), 'AECS') |
                                                            str_detect(toupper(source_name_ch1), 'AIRWAY EPITHELIAL CELLS') ~ 'Airway epithelium',
                                                          
                                                          str_detect(toupper(source_name_ch1), 'ALVEOLAR MACROPHAGE') ~ 'Alveolar macrophage',
                                                          
                                                          str_detect(toupper(source_name_ch1), 'PBMC') |
                                                            str_detect(toupper(source_name_ch1), 'PERIPHERAL BLOOD MONONUCLEAR CELL') ~ 'PBMC',
                                                          
                                                          str_detect(toupper(source_name_ch1), 'BRAIN WHITE MATTER') ~ 'White matter',
                                                          
                                                          str_detect(toupper(source_name_ch1), 'BRONCHIAL') |
                                                            str_detect(toupper(source_name_ch1), 'HBEC') ~ 'Bronchial epithelium',
                                                          
                                                          str_detect(toupper(source_name_ch1), 'BRONCHOALVEOLAR LAVAGE') ~ 'Bronchoalveolar lavage',
                                                          
                                                          str_detect(toupper(source_name_ch1), 'NEUTROPHIL') ~ 'Neutrophil',
                                                          
                                                          str_detect(toupper(source_name_ch1), 'PRIMARY MSCS') ~ 'Mesenchymal stem cell',
                                                          
                                                          str_detect(toupper(source_name_ch1), 'CD14+ MONOCYTES') |
                                                            str_detect(toupper(source_name_ch1), 'THP-1') |
                                                            str_detect(toupper(source_name_ch1), 'MONOCYTES') & 
                                                            !str_detect(toupper(source_name_ch1), 'DENDRITIC') ~ 'Monocyte',
                                                          
                                                          str_detect(toupper(source_name_ch1), 'NK CELLS') ~ 'NK cell',
                                                          
                                                          str_detect(toupper(source_name_ch1), 'ENRICHMENT: CD45+LIN-HLA-DR+') | 
                                                            str_detect(toupper(source_name_ch1), 'MONOCYTE-DERIVED DENDRITIC CELLS') ~ 'Dendritic cell',
                                                          
                                                          str_detect(toupper(source_name_ch1), 'HLADR-CD33+') ~ 'Myeloid-derived suppressor cell',
                                                          
                                                          str_detect(toupper(source_name_ch1), 'THROMBOCYTES') ~ 'Thrombocytes',
                                                          
                                                          str_detect(toupper(source_name_ch1), 'IPSC-DERIVED MOTOR NEURONS') ~ 'IPSC derived neurons',
                                                          
                                                          str_detect(toupper(source_name_ch1), 'PRIMARY NASAL EPITHELIAL CELLS') ~ 'Nasal airway epithelium',
                                                          
                                                          str_detect(toupper(source_name_ch1), 'LUNG SQUAMOUS CELL') ~ 'Lung squamous cell',
                                                          
                                                          str_detect(toupper(source_name_ch1), 'HUH7 CELLS') ~ 'Huh7 cell',
                                                          
                                                          str_detect(toupper(source_name_ch1), 'HL TISSUE') ~ 'Hodgkin Lymphoma',
                                                          
                                                          str_detect(toupper(source_name_ch1), 'PASMC') ~ 'Pulmonary artery smooth muscle',
                                                          
                                                          str_detect(toupper(source_name_ch1), 'NKT') ~ 'Natural killer cell',
                                                          
                                                          str_detect(toupper(source_name_ch1), 'IPSC-DERIVED MOTOR NEURON') ~ 'IPSC-derived neurons',
                                                          
                                                          str_detect(toupper(source_name_ch1), 'NASAL POLYP') ~ 'Nasal polyp',
                                                          
                                                          str_detect(toupper(source_name_ch1), 'BRONCHI') ~ 'Primary bronchial epithelial cells',
                                                          
                                                          #########################################
                                                          
                                                          str_detect(characteristics_ch1, 'cell type: monocyte') |
                                                            str_detect(characteristics_ch1, 'cell type: Monocyte') |
                                                            str_detect(characteristics_ch1, 'cell type: CD14') |
                                                            str_detect(characteristics_ch1, 'Chronic Fatigue Syndrome (FCS),tissue: Adult Human Peripheral Blood,cell type: Monocytes') ~ 'Monocyte',
                                                          
                                                          str_detect(characteristics_ch1, 'cell type: human monocyte-derived macrophages') ~ 'Monocyte-derived macrophage',
                                                          
                                                          str_detect(characteristics_ch1, 'tissue: Brain,lesion type: White matter') ~ 'White matter',
                                                          
                                                          str_detect(characteristics_ch1, 'cell type: Islet of Langerhans') ~ 'Pancreatic islet',
                                                          
                                                          str_detect(characteristics_ch1, 'tissue: airway epithelium') |
                                                            str_detect(characteristics_ch1, 'cell line: A549') |
                                                            str_detect(characteristics_ch1, 'cell type: lung epithelial cells') ~ 'Airway epithelium',
                                                          
                                                          str_detect(characteristics_ch1, ' B cell') ~ 'B cell',
                                                          
                                                          str_detect(characteristics_ch1, 'type: Memory CD4 T cells') |
                                                            str_detect(characteristics_ch1, 'CD3 pos T cells') |
                                                            str_detect(characteristics_ch1, 'cell type: CD3+CD4+CD45RA-CCR7+CD25-CCR4+ T cells') |
                                                              str_detect(characteristics_ch1, 'cell type: T cell') ~ 'T cells',
                                                          
                                                          str_detect(characteristics_ch1, 'tissue: peripheral blood mononuclear cells') | 
                                                            str_detect(characteristics_ch1, 'cell type: PBMC') ~ 'PBMC',
                                                          
                                                          str_detect(characteristics_ch1, 'cell type: iPs derived neurons') |
                                                            str_detect(characteristics_ch1, 'iPSC derived motor neuron') ~ 'IPSC derived neurons',
                                                          
                                                          str_detect(characteristics_ch1, 'sample type: bronchial smooth muscle') ~ 'Airway smooth muscle',
                                                          
                                                          str_detect(characteristics_ch1, 'tissue: tracheobronchial,') ~ 'Bronchia',
                                                          
                                                          str_detect(characteristics_ch1, 'cell type: dermal fibroblasts') |
                                                            str_detect(characteristics_ch1, 'karyotype: Fibroblasts') ~ 'Fibroblast',
                                                          
                                                          str_detect(characteristics_ch1, 'CLL IGHV unmutated') | 
                                                            str_detect(characteristics_ch1, 'unmutated CLL') ~ 'Unmutated CLL',
                                                          
                                                          str_detect(characteristics_ch1, 'CLL IGHV mutated') |
                                                            str_detect(characteristics_ch1, 'Mutated CLL') ~ 'Mutated CLL',
                                                          
                                                          str_detect(characteristics_ch1, 'CD19+ CD3- CD20+ IgD- CD27-') ~ 'CD19+ CD3- CD20+ IgD- CD27-',
                                                          
                                                          str_detect(characteristics_ch1, 'CD19- CD14+ CD3-') ~ 'CD19- CD14+ CD3-',
                                                          
                                                          str_detect(characteristics_ch1, 'fibroblast') ~ 'Fibroblast',
                                                          
                                                          str_detect(characteristics_ch1, 'cell type: neutrophil') ~ 'Neutrophil',
                                                          
                                                          str_detect(characteristics_ch1, 'karyotype: Fibroblasts') ~ 'Fibroblast',
                                                          str_detect(characteristics_ch1, 'karyotype: Fibroblasts') ~ 'Fibroblast',
                                                          str_detect(characteristics_ch1, 'karyotype: Fibroblasts') ~ 'Fibroblast',
                                                          str_detect(characteristics_ch1, 'karyotype: Fibroblasts') ~ 'Fibroblast',
                                                          
                                                          TRUE ~ as.character(.$Tissue_Type)))

#########################################

Selected_Samples <- Selected_Samples %>% mutate(Filtered_Tissue_Type =
                                                  case_when(str_detect(source_name_ch1, 'CD4+') |
                                                              str_detect(Tissue_Type, 'CD4+') |
                                                              str_detect(characteristics_ch1, 'CD4+') ~ 'CD4+ T Cell',
                                                            
                                                            str_detect(source_name_ch1, 'CD8+') |
                                                              str_detect(Tissue_Type, 'CD8+') |
                                                              str_detect(characteristics_ch1, 'CD8+') ~ 'CD8+ T Cell',
                                                            
                                                            str_detect(source_name_ch1, 'TH1') |
                                                              str_detect(Tissue_Type, 'TH1') |
                                                              str_detect(characteristics_ch1, 'TH1') ~ 'Th1 Cell',
                                                            
                                                            str_detect(toupper(source_name_ch1), 'TH2') |
                                                              str_detect(toupper(Tissue_Type), 'TH2') |
                                                              str_detect(toupper(characteristics_ch1), 'TH1') ~ 'Th2 Cell',
                                                            
                                                            str_detect(toupper(source_name_ch1), 'TREG') |
                                                              str_detect(toupper(Tissue_Type), 'TREG') |
                                                              str_detect(toupper(characteristics_ch1), 'TREG') ~ 'Treg Cell',
                                                              
                                                              str_detect(toupper(source_name_ch1), 'MAIT CELLS') ~ 'Mucosal-associated invariant T Cell',
                                                            
                                                            TRUE ~ as.character(.$Filtered_Tissue_Type)))

#########################################

# Remove samples with undefined MeSH terms
Selected_Samples <- Selected_Samples[which(is.na(Selected_Samples$Used_MeSH_Term) != TRUE), ]

setwd("")

save(Selected_Samples, 
     file = 'Selected_Defined_Samples.RData')