import os
import math
import time
import json
import pyreadr
import torch
from torch import nn
import torch.utils.data as td
from sklearn.utils.class_weight import compute_class_weight
import numpy as np

# Check whether cuda is available on device
torch.cuda.is_available()

device = (
    "cuda"
    if torch.cuda.is_available()
    else "mps"
    if torch.backends.mps.is_available()
    else "cpu"
)

# Debug in cpu to get better traceback
# device = "cpu"

# Import data isoalted, aligned, and saved as .RData
Imported_data = pyreadr.read_r('E:/UNC Projects/DNN_Python/Total Dataset/5. MeSH Alignment/Output/Pulmonary_LCPM_Data.RData')

# Print the data keys
print(Imported_data.keys())

# Define the datasets
Training_dataset = Imported_data["Training_Data"]
Testing_dataset = Imported_data["Testing_Data"]

Features = Training_dataset.columns[10:Training_dataset.shape[1]] # lcpm data for each gene
Label = 'Tree_Term_Categorical' # Numeric representation of single tree term

# Map numeric tree term designation to actual tree terms for future export
category_map = dict(zip(Training_dataset['Tree_Term_Categorical'].astype(int), Training_dataset['Tree_Term']))

# Create x (lcpm) and y (tree term categorical) train/test data sets
x_train = Training_dataset[Features].values
y_train = Training_dataset[Label]

x_test = Testing_dataset[Features].values
y_test = Testing_dataset[Label]

print ('Training Set: %d, Test Set: %d \n' % (len(x_train), len(x_test)))

# Set seed
torch.manual_seed(0)

# Create a dataset and loader for the training data and labels
train_x = torch.Tensor(x_train).half().to(device)
train_y = torch.Tensor(y_train).long().to(device)
train_ds = td.TensorDataset(train_x, train_y)
train_loader = td.DataLoader(train_ds, 
                             batch_size = 50,
                             num_workers = 0)

# Create a dataset and loader for the test data and labels
test_x = torch.Tensor(x_test).half().to(device)
test_y = torch.Tensor(y_test).long().to(device)
test_ds = td.TensorDataset(test_x, test_y)
test_loader = td.DataLoader(test_ds, 
                            batch_size = 50,
                            num_workers = 0)

#################################

from DynamicNN import DynamicNN, test, train

layers = int(math.ceil(2/3 * (len(Features) + int(max(Training_dataset[Label]) + 1))))

model = DynamicNN(Input_layer_size = len(Features), 
                  Output_layer_size = int(max(Training_dataset[Label]) + 1), 
                  Number_hidden_layers = 8, 
                  Hidden_layer_neurons = [int(math.ceil(layers/2)),
                                          int(math.ceil(layers/5)), 
                                          int(math.ceil(layers/10)), 
                                          int(math.ceil(layers/10)), 
                                          int(math.ceil(layers/20)), 
                                          int(math.ceil(layers/20)),
                                          int(math.ceil(layers/40)), 
                                          int(math.ceil(layers/80))]).to(device)

print(model)

total_params = sum(p.numel() for p in model.parameters())

total_params

#################################

# Save the model weights
model_file = f'models/Lung_Subset_Model.pt'

model_exists = os.path.exists(model_file)

learning_rate = 0.01
optimizer = torch.optim.SGD(model.parameters(), lr = learning_rate)
optimizer.zero_grad()

if model_exists == True:

    checkpoint = torch.load(model_file, weights_only = True)
    model.load_state_dict(checkpoint['model_state_dict'])
    optimizer.load_state_dict(checkpoint['optimizer_state_dict'])
    epoch_nums = checkpoint['epoch']
    training_loss = checkpoint['loss']

    with open('validation_loss.json', 'r') as file:
        validation_loss = json.load(file)

else:

    epoch_nums = 0
    training_loss = []
    validation_loss = []

#####################################################################

import pandas as pd

All_Classes = pd.concat([Training_dataset[Label]] + [Testing_dataset[Label]])

weights = compute_class_weight(class_weight = 'balanced', 
                               classes = np.unique(y_train), 
                               y = y_train)

train_weights = torch.tensor(weights, dtype = torch.float16).to(device)

# Specify the loss criteria (we'll use CrossEntropyLoss for multi-class classification)
loss_criteria = nn.CrossEntropyLoss(weight = train_weights, reduction = 'mean')

# Train over epochs
epochs = 50
training_epochs = 0

# Calculate the start time
start = time.time()

for epoch in range(1, epochs + 1):

    # print the epoch number
    print('\nEpoch: {}'.format(epoch))

    # Calculate the start time
    train_start = time.time()
    
    # Feed training data into the model to optimize the weights
    train_loss = train(model, 
                       train_loader, 
                       optimizer, 
                       loss_criteria, 
                       autocast_dytpe = torch.float16)
    
    print('\n Training epoch ', epoch, 'took ', round((time.time() - train_start)/60, 2), ' minutes.')

    torch.cuda.empty_cache()

    if epoch % 5 == 0:
        
        # Feed the test data into the model to check its performance
        test_loss = test(model, 
                         test_loader, 
                         loss_criteria, 
                         autocast_dytpe = torch.float16)
        
        validation_loss.append(test_loss)
        training_epochs += 1

    # Log the metrics for this epoch
    epoch_nums += 1
    training_loss.append(train_loss)

# Calculate the end time and time taken
end = time.time()
length = end - start

# Show the results : this can be altered however you like
print("\nTraining over", epoch, "epoch(s) took ", round(length/60, 2), "minutes!")

#################################

%matplotlib inline
from matplotlib import pyplot as plt

plt.plot(list(range(1, epoch_nums + 1)), training_loss)
plt.plot([x * 5 for x in list(range(1, int(epoch_nums/5) + 1))], validation_loss)
plt.xlabel('epoch')
plt.ylabel('loss')
plt.legend(['training', 'validation'], loc='upper right')
plt.show()

#################################

# Set the model to evaluate mode
model.eval()

# Get predictions for the test data
x = torch.Tensor(x_test).float().to(device)
_, predicted = torch.max(model(x).data, 1)

compare = []

for i in range(len(y_test)):
    if y_test[i] == predicted.to("cpu").numpy()[i]:
        compare.append('Correct')
    else:
        compare.append('Incorrect')

Class_key = dict(zip(Testing_dataset['Tree_Term_Categorical'], 
                     Testing_dataset['Tree_Term']))

Actual = []

for i in range(len(y_test)):
    Actual.append(Class_key[y_test[i]])

Predicted = []

for i in range(len(predicted.to("cpu").numpy())):
    Predicted.append(Class_key[predicted.to("cpu").numpy()[i]])

data = {'Actual' : Actual, 
        'Predicted' : Predicted, 
        'Prediction_Result' : compare,
        'Tissue' : Testing_dataset['Tissue_Type'], 
        'Disease_Name' : Testing_dataset['Disease_Name'], 
        'Tree_Term' : Testing_dataset['Tree_Term']}

DNN_predictions = pd.DataFrame(data)
DNN_predictions.to_csv('Output_1_Layer.csv')

#################################

torch.save({
            'epoch': epoch_nums,
            'model_state_dict': model.state_dict(),
            'optimizer_state_dict': optimizer.state_dict(),
            'loss': training_loss,
            }, model_file)

with open('validation_loss.json', 'w') as file:
    json.dump(validation_loss, file)