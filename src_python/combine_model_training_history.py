import glob
import pandas as pd
from operator import itemgetter
import numpy as np

files = glob.glob('../model_training_history/NN_evaluation_seeds1_100_folds5_reducedV3.2_removeN5/*')
print(files)
trainingprogress = {}
validationprogress = {}
trainingprogress_epochs = {}
validationprogress_epochs = {}
indices = [(0, 5), (1, 6), (2, 7), (3, 8), (4, 9)]

for file in files:
    currTrainHistory = pd.read_csv(file, sep='\t', index_col=0)
    for index in indices:
        currTrain = currTrainHistory.iloc[index[0]].dropna()
        currValidation = currTrainHistory.iloc[index[1]].dropna()
        for epoch in range(len(currTrain)):
            percentTrain = round((epoch+1)/len(currTrain), 3)
            if percentTrain not in trainingprogress:
                trainingprogress[percentTrain] = [currTrain[epoch]]
                validationprogress[percentTrain] = [currValidation[epoch]]
            else:
                trainingprogress[percentTrain].append(currTrain[epoch])
                validationprogress[percentTrain].append(currValidation[epoch])
            if epoch not in trainingprogress_epochs:
                trainingprogress_epochs[epoch] = [float(currTrain[epoch])]
                validationprogress_epochs[epoch] = [float(currValidation[epoch])]
            else:
                trainingprogress_epochs[epoch].append(float(currTrain[epoch]))
                validationprogress_epochs[epoch].append(float(currValidation[epoch]))

trainproglist = []
validationproglist = []
trainlist_epochs = []
validationlist_epochs = []

for key, val in trainingprogress.items():
    trainproglist.append((key, np.mean(trainingprogress[key])))

for key, val in validationprogress.items():
    validationproglist.append((key, np.mean(validationprogress[key])))

for key, val in trainingprogress_epochs.items():
    trainlist_epochs.append((key,
                             np.mean(trainingprogress_epochs[key]),
                             np.std(trainingprogress_epochs[key]),
                             len(trainingprogress_epochs[key])))

for key, val in validationprogress_epochs.items():
    validationlist_epochs.append((key,
                                  np.mean(validationprogress_epochs[key]),
                                  np.std(validationprogress_epochs[key]),
                                  len(validationprogress_epochs[key])))

trainproglist = sorted(trainproglist, key=itemgetter(0))
validationproglist = sorted(validationproglist, key=itemgetter(0))
trainlist_epochs = sorted(trainlist_epochs, key=itemgetter(0))
validationlist_epochs = sorted(validationlist_epochs, key=itemgetter(0))

with open('../model_training_history/combinedtraininghistory.txt', 'w+') as f:
    f.write('PercentTrain\ttrainingloss\tvalidationloss\n')
    for i in range(0, len(trainproglist)):
        currentPercent = trainproglist[i][0]
        currentTrainLoss = trainproglist[i][1]
        currentValidationLoss = validationproglist[i][1]
        f.write(str(currentPercent)+'\t'+str(currentTrainLoss)+'\t'+str(currentValidationLoss)+'\n')

with open('../model_training_history/combinedtrainhistory_epochs.txt', 'w+') as f:
    f.write('Epoch\ttrainingloss\tvalidationloss\tstdTrain\tstdValidation\tnum_models\n')
    for i in range(0, len(trainlist_epochs)):
        currentEpoch = trainlist_epochs[i][0]
        currentTrainLoss = trainlist_epochs[i][1]
        currentValidationLoss = validationlist_epochs[i][1]
        stdevTrain = trainlist_epochs[i][2]
        stdevVal = validationlist_epochs[i][2]
        num_Train = trainlist_epochs[i][3]
        f.write(str(currentEpoch) + '\t' +
                str(currentTrainLoss) + '\t' + str(currentValidationLoss) + '\t' +
                str(stdevTrain) + '\t' + str(stdevVal) + '\t' +
                str(num_Train) + '\n')


