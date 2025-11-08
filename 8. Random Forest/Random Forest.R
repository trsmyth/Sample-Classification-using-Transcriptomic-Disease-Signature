rm(list = ls(all.names = TRUE)) # clears global environ

# Load packages
library(tidyverse) # for data cleaning
library(dplyr) 
library(ranger)
library(ConfusionTableR)
library(grid)
library(gridExtra)
library(gridtext)
library(caret)

load('Pulmonary_LCPM_Data.RData')

# Count number of samples in each polarization ID group
group_num <- Training_Data %>% 
  group_by(Tree_Term) %>% 
  summarise(total_count = n(),
            .groups = 'drop')

# Isolate tree term of training data as a factor
Tree_Term <- factor(Training_Data$Tree_Term)

group_num$Tree_Term <- factor(group_num$Tree_Term, 
                              levels = levels(Tree_Term))

# Calculate total number of samples
total <- sum(group_num[, 2])

# Calculate the fraction of whole for each polarization state
fraction <- group_num[, 2]/total

# Calculate ratio for weighing samples
ratio <- 1 - fraction
rownames(ratio) <- group_num$Tree_Term

# Create weight vector with each value set to 0
weight <- as.numeric(Tree_Term)

for (i in 1:length(unique(weight))){
  
  # Replace numeric value representing Tree Term with the ratio for that term
  weight <- replace(weight, weight == i, ratio[i, ])
  
}

Train_meta <- Training_Data[, 1:8]
Training_Data <- Training_Data[, -c(1:6, 8:10)]

Test_meta <- Training_Data[, 1:8]
Testing_Data <- Testing_Data[, -c(1:6, 8:10)]

#############################

# Define mtry for model
best.m = sqrt(ncol(Training_Data) - 1)

Training_Data$Tree_Term <- factor(Training_Data$Tree_Term)

Testing_Data$Tree_Term <- factor(Testing_Data$Tree_Term)

# Run the RF model through ranger package
classifier_RF <- ranger(x = Training_Data[-1], 
                        y = Training_Data$Tree_Term,
                        num.trees = 2000,
                        max.depth = 20, # Limit tree depth to 20 splits if any trees would grow beyond 20 splits
                        mtry = best.m, # mtry set to sqrt(# genes)
                        case.weights = weight, # Weigh model based on sample ratio
                        importance = 'permutation') # Calculate permutation importance

print(classifier_RF)

#######

# Load the group ordering used for the DNN confusion matrix for easier comparison
Group_order <- read.csv("Group_order.csv")

# Reformat group ordering and tree terms to remove punctuation and spacing inconsistencies
Group_order <- gsub(pattern = '[-,+()]', replacement = ' ', Group_order$X0)
Group_order <- gsub(pattern = ' +', replacement = ' ', Group_order)
Testing_Data$Tree_Term <- gsub(pattern = '[-,+()]', replacement = ' ', Testing_Data$Tree_Term)
Testing_Data$Tree_Term <- gsub(pattern = ' +', replacement = ' ', Testing_Data$Tree_Term)

# Predict testing data to determine model performance
p2 <- predict(classifier_RF, Testing_Data) %>% data.frame()
colnames(p2) <- gsub(pattern = '\\.+', replacement = ' ', colnames(p2))

# Reorder prediction results for isolation of top 3 predictions
Top_n_predictions <- t(p2) %>% data.frame() %>%  rownames_to_column()

# Determine the top 3 predictions for each test sample
Top_3_predictions <- lapply(colnames(Top_n_predictions)[2:ncol(Top_n_predictions)], function(Sample){
  
  Single_sample <- order(Top_n_predictions[, Sample], decreasing = TRUE)
  Top_3_predictions <- as.list(Top_n_predictions$rowname[Single_sample][1:3])
  
})

# Create a df with the actual class and each of the top 3 predictions in their own column
Prediction <- data.frame('Actual' = factor(Testing_Data$Tree_Term, levels = Group_order))
Prediction$Prediction = factor(lapply(Top_3_predictions, `[[`, 1), levels = Group_order)
Prediction$Prediction_2 = factor(lapply(Top_3_predictions, `[[`, 2), levels = Group_order)
Prediction$Prediction_3 = factor(lapply(Top_3_predictions, `[[`, 3), levels = Group_order)

#######

# Create a confusion matrix for the actual versus top prediction
confusion_matrix <- as.data.frame(table(Prediction$Actual, Prediction$Prediction))

bin_cm <- ConfusionTableR::binary_class_cm(Prediction$Actual, 
                                           Prediction$Prediction)

# Isolate model performance data
RF_Stats <- data.frame(bin_cm[["confusion_matrix"]][["overall"]][1:5])
colnames(RF_Stats) <- ''

RF_Stats2 <- data.frame(bin_cm[["confusion_matrix"]][["byClass"]])
RF_Stats2 <- RF_Stats2[, 1:7]

# Create tables of performance data
RF_stats_Table <- tableGrob(t(round(RF_Stats, 4))) # Accuracy
RF_stats_Table2 <- tableGrob(round(RF_Stats2, 2)) # Sensitivity/Specificity/Recall

png(file = 'RF_stats_Table.png', height = 100, width = 500)

grid.arrange(RF_stats_Table)

dev.off()

##############

# Plot results of performance data
plot_results <- ggplot(data = confusion_matrix,
                       mapping = aes(x = Var2, y = Var1)) +
  
  geom_tile(aes(fill = Freq),
            show.legend = FALSE) +
  
  geom_text(aes(label = as.character(Freq)), vjust = 1) +
  
  scale_fill_gradient(low = "white",
                      high = "forestgreen", 
                      trans = "log") +
  # Label x and y axis
  labs(title = 'Test Data Confusion Matrix',
       x = "Actual",
       y = "Prediction") + 
  
  # Increase axis text size
  theme(title = element_text(size = 16, face = 'bold'),
        axis.title = element_text(size = 16, face = 'bold'), 
        axis.text = element_text(size = 10), 
        axis.text.x = element_text(angle = 90, hjust = 1))

png(file = 'RF_Confusion_Matrix.png', height = 1000, width = 1000)

print(plot_results)

dev.off()

##############

# Create a confusion matrix for the actual versus top 3 predictions
Top_3_matrix <- data.frame('Actual' = rep('TRUE', nrow(Prediction)), 'Prediction' = NA)

for(i in 1:nrow(Prediction)){
  
  Current_prediction <- as.character(Prediction[i, 'Actual'])
  Top_3_matrix[i, 'Prediction'] <- any(apply(Prediction[i, 2:ncol(Prediction)], 1, function(row) Current_prediction %in% row))
  
}

confusion_matrix <- as.data.frame(table(Top_3_matrix$Actual, Top_3_matrix$Prediction))

# Print the percent of the correct top 3 predictions
print(confusion_matrix[2, 'Freq']/nrow(Top_3_matrix))
      