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
import pandas as pd

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
Imported_data = pyreadr.read_r("Pulmonary_LCPM_Data.RData")

# Print the data keys
print(Imported_data.keys())

# Define the datasets
Training_dataset = Imported_data["Training_Data"]
Testing_dataset = Imported_data["Testing_Data"]

################################################################

Features = Training_dataset.columns[10:Training_dataset.shape[1]] # lcpm data for each gene
Label = 'Tree_Term_Categorical' # Numeric representation of single tree term

# Map numeric tree term designation to actual tree terms for future export
category_map = dict(zip(Training_dataset['Tree_Term_Categorical'].astype(int), Training_dataset['Tree_Term']))

# Create x (lcpm) and y (tree term categorical) train/test data sets
x_train = Training_dataset[Features].values
y_train = Training_dataset[Label]

x_test = Testing_dataset[Features].values
y_test = Testing_dataset[Label]

if np.isnan(x_train).any():
    print("Training data contains NaNs")

if np.isnan(x_test).any():
    print("Testing data contains NaNs")

print('Training Set: %d, Test Set: %d \n' % (len(x_train), len(x_test)))

# Set seed
torch.manual_seed(0)

# Create a dataset and loader for the training data and labels
train_x = torch.Tensor(x_train).half().to(device) # Set count data to half-precision floating-point to save memory
train_y = torch.Tensor(y_train).long().to(device) # Set data label to 64-bit signed integer
train_ds = td.TensorDataset(train_x, train_y)
train_loader = td.DataLoader(train_ds, 
                             batch_size = 100,
                             num_workers = 0)

# Create a dataset and loader for the test data and labels
test_x = torch.Tensor(x_test).half().to(device)
test_y = torch.Tensor(y_test).long().to(device)
test_ds = td.TensorDataset(test_x, test_y)
test_loader = td.DataLoader(test_ds, 
                            batch_size = 100,
                            num_workers = 0)


All_Classes = pd.concat([Training_dataset[Label]] + [Testing_dataset[Label]])

weights = compute_class_weight(class_weight = 'balanced', 
                               classes = np.unique(All_Classes), 
                               y = y_train)

train_weights = torch.tensor(weights, dtype = torch.float16).to(device)

# Specify the loss criteria (we'll use CrossEntropyLoss for multi-class classification)
loss_criteria = nn.CrossEntropyLoss(weight = train_weights, reduction = 'mean')

#################################

from DynamicNN import DynamicNN, test, train

# Define the starting number of neurons per layer. Starting point was determined to be 
# 2/3 the number of features (genes) and labels (disease terms)
layers = int(math.ceil(2/3 * (len(Features) + int(max(Training_dataset[Label]) + 1))))

model = DynamicNN(Input_layer_size = len(Features), 
                  Output_layer_size = int(max(Training_dataset[Label]) + 1), 
                  Number_hidden_layers = 6, 
                  Hidden_layer_neurons = [int(math.ceil(layers)), 
                                          int(math.ceil(layers)/5),
                                          int(math.ceil(layers)/5),
                                          int(math.ceil(layers)/10),
                                          int(math.ceil(layers)/10),
                                          int(math.ceil(layers)/20)]).to(device)

print(model)

total_params = sum(p.numel() for p in model.parameters())

print(total_params)

#################################

# Check if a model checkpoint already exists in case of multiple training runs
# For this analysis, the number of training epochs was tracked and models were
# saved at set numbers of epochs in a folder named '/models'
if os.path.exists(f'models/epoch_nums.npy') == True:
    epoch_nums = np.load(f'models/epoch_nums.npy')

else:
    epoch_nums = 0

# Determine if a model checkpoint is present
model_file = f'models/Lung_Subset_Model_{epoch_nums}.pt'

model_exists = os.path.exists(model_file)

# If a checkpoint does exist, load the corresponding data from that checkpoint
if model_exists == True:

    training_epochs = np.load(f'models/training_epochs_{epoch_nums}.npy')
    percent_correct = np.load(f'models/percent_correct_{epoch_nums}.npy').tolist()
    percent_top_3_correct = np.load(f'models/percent_top_3_correct_{epoch_nums}.npy').tolist()
    top_3_predicted = np.load(f'models/top_3_predicted_{epoch_nums}.npy').tolist()
    
    learning_rate = 0.005
    optimizer = torch.optim.SGD(model.parameters(), lr = learning_rate)
    optimizer.zero_grad()

    checkpoint = torch.load(model_file, weights_only = False)
    model.load_state_dict(checkpoint['model_state_dict'])
    optimizer.load_state_dict(checkpoint['optimizer_state_dict'])
    epoch_nums = checkpoint['epoch']
    training_loss = checkpoint['loss']

    with open(f'models/validation_loss_{epoch_nums}.json', 'r') as file:
        validation_loss = json.load(file)

    del checkpoint 
    torch.cuda.empty_cache()

else:

    training_epochs = 0
    training_loss = []
    validation_loss = []
    percent_correct = []
    percent_top_3_correct = [] 
    top_3_predicted = []

    learning_rate = 0.005
    optimizer = torch.optim.SGD(model.parameters(), lr = learning_rate)
    optimizer.zero_grad()

#####################################################################

# Train over epochs
epochs = 500

# Calculate the start time
start = time.time()

for epoch in range(1, epochs + 1):

    # Add 1 to epoch counter
    epoch_nums += 1

    # print the epoch number
    print('\nEpoch: {}'.format(epoch_nums))

    # Calculate the start time
    train_start = time.time()
    
    # Feed training data into the model to optimize the weights
    train_loss = train(model, 
                       train_loader, 
                       optimizer, 
                       loss_criteria, 
                       autocast_dytpe = torch.float16)
    
    print('\n Training epoch', epoch, 'took', round((time.time() - train_start)/60, 2), 'minutes.')

    torch.cuda.empty_cache()

    if epoch_nums % 5 == 0:
        
        # Feed the test data into the model to check its performance
        test_loss = test(model, 
                         test_loader, 
                         loss_criteria, 
                         autocast_dytpe = torch.float16, 
                         top_prediction_cutoff = 3)

        validation_loss.append(test_loss[0])
        percent_correct.append(test_loss[1])
        percent_top_3_correct.append(test_loss[2])
        top_3_predicted.append(test_loss[3])
        training_epochs += 1

    # Log the metrics for this epoch
    training_loss.append(train_loss)

    # Save the model every 250 epochs
    if epoch_nums % 250 == 0:

        np.save(f'models/epoch_nums', epoch_nums)
        np.save(f'models/training_epochs_{epoch_nums}', training_epochs)
        np.save(f'models/Percent_Top_3_Correct_{epoch_nums}', percent_top_3_correct)
        np.save(f'models/percent_correct_{epoch_nums}', percent_correct)
        np.save(f'models/top_3_predicted_{epoch_nums}', top_3_predicted)

        # Save the model weights
        model_file = f'models/Lung_Subset_Model_{epoch_nums}.pt'

        torch.save({
                    'epoch': epoch_nums,
                    'model_state_dict': model.state_dict(),
                    'optimizer_state_dict': optimizer.state_dict(),
                    'loss': training_loss,
                    }, model_file)

        with open(f'models/validation_loss_{epoch_nums}.json', 'w') as file:
            json.dump(validation_loss, file)

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
plt.legend(['training', 'validation'], loc = 'upper right')
plt.show()

plt.plot([x * 5 for x in list(range(1, int(epoch_nums/5) + 1))], percent_correct)
plt.plot([x * 5 for x in list(range(1, int(epoch_nums/5) + 1))], percent_top_3_correct)
plt.xlabel('epoch')
plt.ylabel('Percent Correct')
plt.ylim(0, 100)
plt.legend(['Percent Correct', 'Percent Top 3 Correct'], loc = 'upper right')
plt.show()

#################################

Control_number = len(np.unique(Testing_dataset['Tree_Term'][Testing_dataset['Tree_Term'].str.contains('CONTROL')].astype('category').cat.codes))
Total_number = len(np.unique(Testing_dataset['Tree_Term'].astype('category').cat.codes))
array_map = np.vectorize(lambda x: category_map.get(x, x))

from sklearn.metrics import confusion_matrix

# Set the model to evaluate mode
model.eval()

# Get predictions for the test data
x = torch.Tensor(test_x).float().to(device)
_, predicted = torch.max(model(x).data, 1)

predicted = predicted.to("cpu").numpy()
test_y = test_y.to('cpu')

# Plot the confusion matrix
cm = confusion_matrix(predicted, test_y)
cm = np.log10(cm)

plt.figure(figsize = (10, 10))
plt.imshow(cm[::-1], interpolation = "nearest", cmap = plt.cm.Blues)
plt.colorbar()

tick_marks = np.arange(int(max(Testing_dataset[Label].values) + 1), step = 1)
text_marks = array_map(tick_marks)

pd.DataFrame(text_marks).to_csv('Group_order.csv')

plt.yticks(tick_marks[::-1], text_marks)
plt.xticks(tick_marks, text_marks, rotation = 90)

plt.axline((0, int(max(Testing_dataset[Label].values))), 
           slope = -1, 
           color = 'blue', 
           linestyle = '--')

plt.axvline(x = Control_number - 0.5, 
            ymin = 0, 
            ymax = Control_number/Total_number, 
            color = 'red', 
            linestyle = '--', 
            linewidth = 1)

plt.axhline(y = Total_number - Control_number - 0.5, 
            xmin = 0, 
            xmax = Control_number/Total_number, 
            color = 'red', 
            linestyle = '--', 
            linewidth = 1)

plt.ylabel("Predicted Species")
plt.xlabel("Actual Species")
plt.show()

cm_prediction_file = f'models/cm_{epoch_nums}.csv'

cm_predictions = pd.DataFrame(cm)
cm_predictions.to_csv(cm_prediction_file)

############

Softmax_appliaction = torch.nn.Softmax(dim = 1)
output = Softmax_appliaction(model(x))
output = output.to('cpu')
top_prediction_number = 3

def get_top_x_predictions(tensor, x):

  top_x_values, top_x_indices = torch.topk(tensor, x)

  return top_x_values, top_x_indices

top_percentages = get_top_x_predictions(output, top_prediction_number)

####

compare = []

for i in range(len(test_y)):
    if test_y[i] == predicted[i]:
        compare.append('Correct')
    else:
        compare.append('Incorrect')

####

Top_3 = []

for i in range(len(test_y)):
    if test_y[i] in top_percentages[1][i]:
        Top_3.append('Yes')
    else:
        Top_3.append('No')

####

Class_key = dict(zip(Testing_dataset['Tree_Term_Categorical'], 
                     Testing_dataset['Tree_Term']))

####

Actual = []

for i in range(len(test_y)):
    Actual.append(Class_key[int(test_y[i])])

####

Predicted = []

for i in range(len(predicted)):
    Predicted.append(Class_key[predicted[i]])

####

top_predictions = []

for i in range(len(top_percentages[0])):
    
    current_value = 100 * top_percentages[0][i].detach().numpy()
    current_key = top_percentages[1][i].numpy()

    Combined_data = []

    for i in range(len(current_key)):

        Combined_data.append((round(current_value[i].item(), 2), Class_key[current_key[i]]))

    top_predictions.append(Combined_data)

####

data = {'geo_accession' : Testing_dataset['geo_accession'], 
        'Prediction_Result' : compare,
        'Actual_in_Top_3' : Top_3,
        'Actual' : Actual, 
        'Top_Predictions' : top_predictions, 
        'Disease_Name' : Testing_dataset['Disease_Name']}

DNN_prediction_file = f'models/DNN_output_epoch_{epoch_nums}.csv'

DNN_predictions = pd.DataFrame(data)
DNN_predictions.to_csv(DNN_prediction_file)