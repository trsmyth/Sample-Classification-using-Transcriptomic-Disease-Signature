import torch
from torch import nn

class DynamicNN(nn.Module):
    def __init__(self, Input_layer_size, Output_layer_size, Number_hidden_layers, Hidden_layer_neurons):
        super(DynamicNN, self).__init__()
        self.Input_layer_size = Input_layer_size

        if type(self.Input_layer_size) != int:
            raise Exception("Input_layer_size reqires a single integer representing the number of input neurons")

        self.Output_layer_size = Output_layer_size

        if type(self.Output_layer_size) != int:
            raise Exception("Output_layer_size reqires a single integer representing the number of output neurons")
        
        self.Number_hidden_layers = Number_hidden_layers
        self.Hidden_layer_neurons = Hidden_layer_neurons

        if len(self.Hidden_layer_neurons) != self.Number_hidden_layers:
            raise Exception("Length of Hidden_layer_neurons must match Number_hidden_layers")

        self.Input_layer = nn.Linear(self.Input_layer_size, self.Hidden_layer_neurons[0])
        self.Hidden_layers = nn.ModuleList()

        if self.Number_hidden_layers == 1:

            # Create the output layer
            self.Output_layer = nn.Linear(self.Input_layer.out_features, self.Output_layer_size)

        else:

            # Create hidden layers
            for i in range(1, self.Number_hidden_layers):
                
                if i == 1:
                    self.Hidden_layers.append(nn.Linear(self.Input_layer.out_features, self.Hidden_layer_neurons[i]))
                    
                else:
                    # Create hidden layer with input equal to previous layer output and 
                    # output equal to the next hidden_layer_neurons number
                    self.Hidden_layers.append(nn.Linear(self.Hidden_layers[i - 2].out_features, self.Hidden_layer_neurons[i]))
                    
            # Create the output layer
            self.Output_layer = nn.Linear(self.Hidden_layers[len(self.Hidden_layers) - 1].out_features, self.Output_layer_size)

    def forward(self, x):
        x = torch.relu(self.Input_layer(x))

        for layer in self.Hidden_layers:
            x = torch.relu(layer(x))

        x = self.Output_layer(x)

        return x

####################################################

def train(model, data_loader, optimizer, loss_criteria, autocast_dytpe = torch.float16):

    # Set the model to training mode
    model.train()
    train_loss = 0
    
    for batch, tensor in enumerate(data_loader):

        data, target = tensor

        # feed forward
        optimizer.zero_grad()
        
        with torch.autocast(device_type = 'cuda', dtype = autocast_dytpe):
            out = model(data)
            loss = loss_criteria(out, target)
            train_loss += loss.item()
            
        # backpropagate
        loss.backward()

        nn.utils.clip_grad_norm_(model.parameters(), 1.0)

        optimizer.step()
        del loss, out

    # Return average loss
    avg_loss = train_loss / (batch + 1)
    print('\nTraining set: Average loss: {:.6f}'.format(avg_loss))
    return avg_loss

#################################

import numpy as np
from functorch.experimental.control_flow import map

def test(model, data_loader, loss_criteria, autocast_dytpe = torch.float16, top_prediction_cutoff = 3):
    
    # Switch the model to evaluation mode (so we don't backpropagate)
    model.eval()
    test_loss = 0
    correct = 0
    top_3_correct = 0
    top_3_predicted = []

    ##########################

    Softmax_appliaction = torch.nn.Softmax(dim = 1)

    def get_top_predictions(tensor, x):

        """With an input of a tensor and an integer x, find the top x softmax
        values and the index of the top x softmax values which corresponds to
        the top x predicted classes"""
        top_x_values, top_x_indices = torch.topk(tensor, x)

        return top_x_values, top_x_indices
    
    ##########################

    with torch.no_grad():

        batch_count = 0

        for batch, tensor in enumerate(data_loader):
            batch_count += 1
            data, target = tensor
            
            with torch.autocast(device_type = 'cuda', dtype = autocast_dytpe):

                # Get the predictions
                out = model(data)

            # calculate the loss
            test_loss += loss_criteria(out, target).item()

            # Calculate the accuracy
            _, predicted = torch.max(out.data, 1)
            correct += torch.sum(target == predicted).item()
            
            output = Softmax_appliaction(out)

            top_predictions = get_top_predictions(output, top_prediction_cutoff)

            for i in range(len(output)):

                top_3_predicted.append(top_predictions[1][i].tolist())

                if target[i] in top_predictions[1][i]:
                    top_3_correct += 1
                    
            # Calculate the top 3 prediction accuracy
            del out

    ##########################
            
    # Calculate the average loss and total accuracy for this epoch
    avg_loss = test_loss / batch_count
    percent_correct = 100. * correct / len(data_loader.dataset)
    percent_top_3_correct = 100. * top_3_correct / len(data_loader.dataset)

    print('\nValidation set: Average loss: {:.6f}, Accuracy: {}/{} ({:.0f}%)\n'.format(
    avg_loss, correct, len(data_loader.dataset), percent_correct))
        
    ##########################
    
    # return average loss for the epoch
    return avg_loss, percent_correct, percent_top_3_correct, top_3_predicted