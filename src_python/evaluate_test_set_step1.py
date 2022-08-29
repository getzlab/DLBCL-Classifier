import re
import pandas as pd
import numpy as np
import os
import random
import classify_generic as CG
import correlation_plot as CP
import sklearn.metrics
import calculate_model_metrics as CM
import format_data as FD

DEFAULT_SEED = 123
NUM_CORRELATION_ITERS = 1000
WINDOW_SIZE = 10
STEP_SIZE = 10

random.seed(DEFAULT_SEED)
np.random.seed(DEFAULT_SEED)

datafile = "../data_tables/gsm/DLBCL_Staudt_Shipp_CL.for_classifier_training.classifier_subset.fix_sv.fix_ploidy.17-Aug-2022.txt"
targetfile = "../data_tables/confidence_tables/baseline_probabilities.connectivity_based.sensitivity_power2.Aug_17_2022.tsv"
qval_file = '../data_tables/qval_dfs/fisher_exact_5x2_17-Aug-2022.combined.tsv'
qval_df = pd.read_csv(qval_file, sep='\t', index_col=0)
sig_genes = list(qval_df.loc[qval_df['q'] < 0.10].index)

testSamples = list(pd.read_csv("../data_tables/train_test_sets/TestingSet_149Subset_May2021.txt",
                                    sep='\t', header=None, index_col=0).index)


labels = pd.read_csv(targetfile, sep='\t', index_col=0)

data = pd.read_csv(datafile, sep='\t', index_col=0)
data = data[testSamples].T
data = data[sig_genes]

reduced_data = FD.construct_reduced_winning_version(data)
reduced_data.to_csv('../evaluation_test_set/test_set_reducedGSM.tsv', sep='\t')

#labels = pd.DataFrame(labels.idxmax(axis=1).map({'C1': 1, 'C2': 2, 'C3': 3, 'C4': 4, 'C5': 5}))
#labels.columns = ['cluster']

#########################
# Classify the test set #
#########################
predDF_unsorted = CG.classify_samples_winning_model(reduced_data)
labels_test = labels.loc[predDF_unsorted.index]

predDF_unsorted['TrueCluster'] = labels_test['cluster']
predDF_unsorted['Correctness'] = np.equal(predDF_unsorted['PredictedCluster'], predDF_unsorted['TrueCluster'])

predDF = predDF_unsorted.sort_values(by='Confidence')

###########################
# Confusion matrices calc #
###########################

top70DF = predDF.iloc[int(predDF.shape[0]*0.30)::]
above90DF = predDF.loc[predDF['Confidence'] >= 0.90]
b_80_90DF = predDF.loc[(predDF['Confidence'] >= 0.80) & (predDF['Confidence'] < 0.90)]
b_70_80DF = predDF.loc[(predDF['Confidence'] > 0.70) & (predDF['Confidence'] < 0.80)]
below70DF = predDF.loc[predDF['Confidence'] <= 0.70]

confusionMatrixTop = sklearn.metrics.confusion_matrix(top70DF['TrueCluster'], top70DF['PredictedCluster'])
accuracyTop70 = np.trace(confusionMatrixTop)
accuracyTop70 = accuracyTop70 / np.sum(confusionMatrixTop)

confusionMatrixAll = sklearn.metrics.confusion_matrix(predDF['TrueCluster'], predDF['PredictedCluster'])
allAccuracy = np.trace(confusionMatrixAll)
allAccuracy = allAccuracy / np.sum(confusionMatrixAll)

confusionMatrixAbove90 = sklearn.metrics.confusion_matrix(above90DF['TrueCluster'], above90DF['PredictedCluster'])
accuracyAbove90 = np.trace(confusionMatrixAbove90) / np.sum(confusionMatrixAbove90)

confusionMatrix_b_80_90 = sklearn.metrics.confusion_matrix(b_80_90DF['TrueCluster'], b_80_90DF['PredictedCluster'])
accuracy_b_80_90 = np.trace(confusionMatrix_b_80_90) / np.sum(confusionMatrix_b_80_90)

confusionMatrix_b_70_80 = sklearn.metrics.confusion_matrix(b_70_80DF['TrueCluster'], b_70_80DF['PredictedCluster'])
accuracy_b_70_80 = np.trace(confusionMatrix_b_70_80) / np.sum(confusionMatrix_b_70_80)

confusionMatrixBelow70 = sklearn.metrics.confusion_matrix(below70DF['TrueCluster'], below70DF['PredictedCluster'])
accuracyBelow70 = np.trace(confusionMatrixBelow70) / np.sum(confusionMatrixBelow70)


print(confusionMatrixAll)
print('Accuracy:', allAccuracy)
print(confusionMatrixTop)
print('Accuracy Top 70%:', accuracyTop70)
print(confusionMatrixAbove90)
print('Accuracy Above 0.90 Confidence:', accuracyAbove90)
print(confusionMatrix_b_80_90)
print('Accuracy Between 0.80, 0.90:', accuracy_b_80_90)
print(confusionMatrix_b_70_80)
print('Accuracy Between 0.70, 0.80:', accuracy_b_70_80)
print(confusionMatrixBelow70)
print('Accuracy Below 0.70 Confidence:', accuracyBelow70)

with open('../data_tables/confusion_matrices/test_set/confusionmatrices' + '.txt', 'w+') as confF:
    confF.write('Format: TRUE labels in ROWS, PREDICTED labels in COLUMNS\n')
    confF.write('Confusion Matrix All\n')
    writestr = np.array2string(confusionMatrixAll, separator='\t')
    writestr = writestr.replace('[', ' ')
    writestr = writestr.replace(']', ' ')
    confF.write(writestr)
    confF.write('\nOverall Accuracy: ' + str(allAccuracy))
    confF.write('\n')
    confF.write('Confusion Matrix Top 70%\n')
    writestr = np.array2string(confusionMatrixTop, separator='\t')
    writestr = writestr.replace('[', ' ')
    writestr = writestr.replace(']', ' ')
    confF.write(writestr)
    confF.write('\nAccuracy Top 70%: ' + str(accuracyTop70))
    confF.write('\n')
    confF.write('Confusion Matrix Above 0.90\n')
    writestr = np.array2string(confusionMatrixAbove90, separator='\t')
    writestr = writestr.replace('[', ' ')
    writestr = writestr.replace(']', ' ')
    confF.write(writestr)
    confF.write('\nAccuracy Above 0.90: ' + str(accuracyAbove90))
    confF.write('\n')
    confF.write('Confusion Matrix 0.80_0.90\n')
    writestr = np.array2string(confusionMatrix_b_80_90, separator='\t')
    writestr = writestr.replace('[', ' ')
    writestr = writestr.replace(']', ' ')
    confF.write(writestr)
    confF.write('\nAccuracy 0.80 < Conf <= 0.90: ' + str(accuracy_b_80_90))
    confF.write('\n')
    confF.write('Confusion Matrix 0.70_0.80\n')
    writestr = np.array2string(confusionMatrix_b_70_80, separator='\t')
    writestr = writestr.replace('[', ' ')
    writestr = writestr.replace(']', ' ')
    confF.write(writestr)
    confF.write('\nAccuracy 0.70 < Conf <= 0.80: ' + str(accuracy_b_70_80))
    confF.write('\n')
    confF.write('Confusion Matrix <= 0.70\n')
    writestr = np.array2string(confusionMatrixBelow70, separator='\t')
    writestr = writestr.replace('[', ' ')
    writestr = writestr.replace(']', ' ')
    confF.write(writestr)
    confF.write('\nAccuracy Conf <= 0.70: ' + str(accuracyBelow70))

lrx = np.array(predDF['Confidence'])
lry = np.array(predDF['Correctness'])

#####################
# Metric evaluation #
#####################

allAccuracy, lowerAcc, upperAcc = CM.modelAccuracy(predDF_unsorted['PredictedCluster'], predDF_unsorted['TrueCluster'])
Kappa, lowerKappa, upperKappa = CM.modelKappa(predDF_unsorted['Confidence'], np.equal(predDF_unsorted['PredictedCluster'], predDF_unsorted['TrueCluster']))
performance, lowerPerf, upperPerf = CM.modelPerformance(allAccuracy, Kappa, lowerAcc, upperAcc, lowerKappa, upperKappa)

# create plots w/ no bootstraps
_, _, _, _ = CP.xyresiduals_window_weighted_fractionPerWindow(  lrx, lry, 'None', seed=123, bootstrapwindows=False,
                                                                 windowsize=WINDOW_SIZE, step=STEP_SIZE, computeplots=True,
                                                                 saveFileName='../plots/test_set/NN_reducedV3.2_removeN5_nfeatures21_testset',
                                                                 format='png')
_, _, _, _ = CP.xyresiduals_window_weighted_fractionPerWindow(  lrx, lry, 'None', seed=123, bootstrapwindows=False,
                                                                windowsize=WINDOW_SIZE, step=STEP_SIZE, computeplots=True,
                                                                saveFileName='../plots/test_set/NN_reducedV3.2_removeN5_nfeatures21_testset',
                                                                format='pdf')

print('Final Metrics\n Accuracy:', allAccuracy, '\n',
      'Kappa:', Kappa, '\n',
      'Performance:', performance, '\n',
      'Mean confidence:', np.mean(predDF['Confidence']), '\n',
      '1-sigma confidence:', np.std(predDF['Confidence']), '\n',)

base_fn = 'NN_reducedV3.2_removeN5_'
app_fn = 'nfeatures21_testsetEval'
ext = '.tsv'

model_evals = {base_fn + app_fn: {'Accuracy': allAccuracy, 'Kappa': Kappa, 'Performance': performance,
                                  'lowerAcc': lowerAcc, 'upperAcc': upperAcc,
                                  'lowerKappa': lowerKappa, 'upperKappa': upperKappa,
                                  'lowerPerf': lowerPerf, 'upperPerf': upperPerf}}

predDF.to_csv("../evaluation_test_set/" + base_fn + app_fn + ext, sep='\t')

# set up reduced frames
noarms_reduced, _ = FD.format_inputs(datafile, targetfile, testSamples, reduced_version='3.2', remove_largest_n=5,
                                     no_arms=True, drop_empty_vectors=False)
nofocals_reduced, _ = FD.format_inputs(datafile, targetfile, testSamples, reduced_version='3.2', remove_largest_n=5,
                                       no_focals=True, drop_empty_vectors=False)
noscnas_reduced, _ = FD.format_inputs(datafile, targetfile, testSamples, reduced_version='3.2', remove_largest_n=5,
                                      no_cna=True, drop_empty_vectors=False)
nosv_reduced, _ = FD.format_inputs(datafile, targetfile, testSamples, reduced_version='3.2', remove_largest_n=5,
                                   no_sv=True, drop_empty_vectors=False)
nomut_reduced, _ = FD.format_inputs(datafile, targetfile, testSamples, reduced_version='3.2', remove_largest_n=5,
                                    no_muts=True, drop_empty_vectors=False)
ploidy_reduced, _ = FD.format_inputs(datafile, targetfile, testSamples, reduced_version='3.2', remove_largest_n=5,
                                     ploidy=True, drop_empty_vectors=False)
coo_reduced, _ = FD.format_inputs(datafile, targetfile, testSamples, reduced_version='3.2', remove_largest_n=5,
                                  coo=True, drop_empty_vectors=False)
coo_ploidy_reduced, _ = FD.format_inputs(datafile, targetfile, testSamples, reduced_version='3.2', remove_largest_n=5,
                                         coo=True, ploidy=True, drop_empty_vectors=False)
# no_bcl6_reduced, _ = FD.format_inputs(datafile, targetfile, testSamples, reduced_version='3.1', qval=0.05, remove_largest_n=5,
#                                       nosvbcl6=True, drop_empty_vectors=False)


models = [noarms_reduced, nofocals_reduced, noscnas_reduced, nosv_reduced, nomut_reduced]
names = ['no.arms', 'no.focals', 'no.scnas', 'no.sv', 'no.muts']

for i in range(len(models)):
    mdl = models[i]
    name = names[i]

    predDF_unsorted = CG.classify_samples_winning_model(mdl)
    labels_test = labels.loc[predDF_unsorted.index]

    predDF_unsorted['TrueCluster'] = labels_test['cluster']
    predDF_unsorted['Correctness'] = np.equal(predDF_unsorted['PredictedCluster'], predDF_unsorted['TrueCluster'])

    predDF = predDF_unsorted.sort_values(by='Confidence')

    top70DF = predDF.iloc[int(predDF.shape[0] * 0.30)::]
    above90DF = predDF.loc[predDF['Confidence'] >= 0.90]
    b_80_90DF = predDF.loc[(predDF['Confidence'] >= 0.80) & (predDF['Confidence'] < 0.90)]
    b_70_80DF = predDF.loc[(predDF['Confidence'] > 0.70) & (predDF['Confidence'] < 0.80)]
    below70DF = predDF.loc[predDF['Confidence'] <= 0.70]

    confusionMatrixTop = sklearn.metrics.confusion_matrix(top70DF['TrueCluster'], top70DF['PredictedCluster'])
    accuracyTop70 = np.trace(confusionMatrixTop)
    accuracyTop70 = accuracyTop70 / np.sum(confusionMatrixTop)

    confusionMatrixAll = sklearn.metrics.confusion_matrix(predDF['TrueCluster'], predDF['PredictedCluster'])
    allAccuracy = np.trace(confusionMatrixAll)
    allAccuracy = allAccuracy / np.sum(confusionMatrixAll)

    confusionMatrixAbove90 = sklearn.metrics.confusion_matrix(above90DF['TrueCluster'], above90DF['PredictedCluster'])
    accuracyAbove90 = np.trace(confusionMatrixAbove90) / np.sum(confusionMatrixAbove90)

    confusionMatrix_b_80_90 = sklearn.metrics.confusion_matrix(b_80_90DF['TrueCluster'], b_80_90DF['PredictedCluster'])
    accuracy_b_80_90 = np.trace(confusionMatrix_b_80_90) / np.sum(confusionMatrix_b_80_90)

    confusionMatrix_b_70_80 = sklearn.metrics.confusion_matrix(b_70_80DF['TrueCluster'], b_70_80DF['PredictedCluster'])
    accuracy_b_70_80 = np.trace(confusionMatrix_b_70_80) / np.sum(confusionMatrix_b_70_80)

    confusionMatrixBelow70 = sklearn.metrics.confusion_matrix(below70DF['TrueCluster'], below70DF['PredictedCluster'])
    accuracyBelow70 = np.trace(confusionMatrixBelow70) / np.sum(confusionMatrixBelow70)

    with open('../data_tables/confusion_matrices/test_set/confusionmatrices_' + name + '.txt', 'w+') as confF:
        confF.write('Format: TRUE labels in ROWS, PREDICTED labels in COLUMNS\n')
        confF.write('Confusion Matrix All\n')
        writestr = np.array2string(confusionMatrixAll, separator='\t')
        writestr = writestr.replace('[', ' ')
        writestr = writestr.replace(']', ' ')
        confF.write(writestr)
        confF.write('\nOverall Accuracy: ' + str(allAccuracy))
        confF.write('\n')
        confF.write('Confusion Matrix Top 70%\n')
        writestr = np.array2string(confusionMatrixTop, separator='\t')
        writestr = writestr.replace('[', ' ')
        writestr = writestr.replace(']', ' ')
        confF.write(writestr)
        confF.write('\nAccuracy Top 70%: ' + str(accuracyTop70))
        confF.write('\n')
        confF.write('Confusion Matrix Above 0.90\n')
        writestr = np.array2string(confusionMatrixAbove90, separator='\t')
        writestr = writestr.replace('[', ' ')
        writestr = writestr.replace(']', ' ')
        confF.write(writestr)
        confF.write('\nAccuracy Above 0.90: ' + str(accuracyAbove90))
        confF.write('\n')
        confF.write('Confusion Matrix 0.80_0.90\n')
        writestr = np.array2string(confusionMatrix_b_80_90, separator='\t')
        writestr = writestr.replace('[', ' ')
        writestr = writestr.replace(']', ' ')
        confF.write(writestr)
        confF.write('\nAccuracy 0.80 < Conf <= 0.90: ' + str(accuracy_b_80_90))
        confF.write('\n')
        confF.write('Confusion Matrix 0.70_0.80\n')
        writestr = np.array2string(confusionMatrix_b_70_80, separator='\t')
        writestr = writestr.replace('[', ' ')
        writestr = writestr.replace(']', ' ')
        confF.write(writestr)
        confF.write('\nAccuracy 0.70 < Conf <= 0.80: ' + str(accuracy_b_70_80))
        confF.write('\n')
        confF.write('Confusion Matrix <= 0.70\n')
        writestr = np.array2string(confusionMatrixBelow70, separator='\t')
        writestr = writestr.replace('[', ' ')
        writestr = writestr.replace(']', ' ')
        confF.write(writestr)
        confF.write('\nAccuracy Conf <= 0.70: ' + str(accuracyBelow70))

    allAccuracy, lowerAcc, upperAcc = CM.modelAccuracy(predDF_unsorted['PredictedCluster'], predDF_unsorted['TrueCluster'])
    Kappa, lowerKappa, upperKappa = CM.modelKappa(predDF_unsorted['Confidence'], np.equal(predDF_unsorted['PredictedCluster'], predDF_unsorted['TrueCluster']))
    performance, lowerPerf, upperPerf = CM.modelPerformance(allAccuracy, Kappa, lowerAcc, upperAcc, lowerKappa, upperKappa)

    model_name = base_fn + name + '_' + app_fn
    print(model_name)
    model_evals[model_name] = {'Accuracy': allAccuracy, 'Kappa': Kappa, 'Performance': performance,
                               'lowerAcc': lowerAcc, 'upperAcc': upperAcc,
                               'lowerKappa': lowerKappa, 'upperKappa': upperKappa,
                               'lowerPerf': lowerPerf, 'upperPerf': upperPerf}

    predDF.to_csv("../evaluation_test_set/preds_" + model_name + ext, sep='\t')

# Extra classifications for other trained models with # features different (+coo, +ploidy)
frames = [ploidy_reduced, coo_reduced, coo_ploidy_reduced]
models_2 = ['NN_evaluation_seeds1_100_folds5_reducedV3.2_coo_removeN5',
            'NN_evaluation_seeds1_100_folds5_ploidy_reducedV3.2_removeN5',
            'NN_evaluation_seeds1_100_folds5_ploidy_reducedV3.2_coo_removeN5']
model_names = ['NN_reducedV3.2_removeN5_coo_nfeatures22_testsetEval',
               'NN_reducedV3.2_removeN5_ploidy_nfeatures22_testsetEval',
               'NN_reducedV3.2_removeN5_coo.ploidy_nfeatures23_testsetEval']

for i in range(len(models_2)):
    fr = frames[i]
    mdl = models_2[i]
    name = model_names[i]

    predDF_unsorted = CG.classify_samples_generic(fr, mdl)
    labels_test = labels.loc[predDF_unsorted.index]

    predDF_unsorted['TrueCluster'] = labels_test['cluster']
    predDF_unsorted['Correctness'] = np.equal(predDF_unsorted['PredictedCluster'], predDF_unsorted['TrueCluster'])

    predDF = predDF_unsorted.sort_values(by='Confidence')

    top70DF = predDF.iloc[int(predDF.shape[0] * 0.30)::]
    above90DF = predDF.loc[predDF['Confidence'] >= 0.90]
    b_80_90DF = predDF.loc[(predDF['Confidence'] >= 0.80) & (predDF['Confidence'] < 0.90)]
    b_70_80DF = predDF.loc[(predDF['Confidence'] > 0.70) & (predDF['Confidence'] < 0.80)]
    below70DF = predDF.loc[predDF['Confidence'] <= 0.70]

    confusionMatrixTop = sklearn.metrics.confusion_matrix(top70DF['TrueCluster'], top70DF['PredictedCluster'])
    accuracyTop70 = np.trace(confusionMatrixTop)
    accuracyTop70 = accuracyTop70 / np.sum(confusionMatrixTop)

    confusionMatrixAll = sklearn.metrics.confusion_matrix(predDF['TrueCluster'], predDF['PredictedCluster'])
    allAccuracy = np.trace(confusionMatrixAll)
    allAccuracy = allAccuracy / np.sum(confusionMatrixAll)

    confusionMatrixAbove90 = sklearn.metrics.confusion_matrix(above90DF['TrueCluster'], above90DF['PredictedCluster'])
    accuracyAbove90 = np.trace(confusionMatrixAbove90) / np.sum(confusionMatrixAbove90)

    confusionMatrix_b_80_90 = sklearn.metrics.confusion_matrix(b_80_90DF['TrueCluster'], b_80_90DF['PredictedCluster'])
    accuracy_b_80_90 = np.trace(confusionMatrix_b_80_90) / np.sum(confusionMatrix_b_80_90)

    confusionMatrix_b_70_80 = sklearn.metrics.confusion_matrix(b_70_80DF['TrueCluster'], b_70_80DF['PredictedCluster'])
    accuracy_b_70_80 = np.trace(confusionMatrix_b_70_80) / np.sum(confusionMatrix_b_70_80)

    confusionMatrixBelow70 = sklearn.metrics.confusion_matrix(below70DF['TrueCluster'], below70DF['PredictedCluster'])
    accuracyBelow70 = np.trace(confusionMatrixBelow70) / np.sum(confusionMatrixBelow70)

    with open('../data_tables/confusion_matrices/test_set/confusionmatrices_' + name + '.txt', 'w+') as confF:
        confF.write('Format: TRUE labels in ROWS, PREDICTED labels in COLUMNS\n')
        confF.write('Confusion Matrix All\n')
        writestr = np.array2string(confusionMatrixAll, separator='\t')
        writestr = writestr.replace('[', ' ')
        writestr = writestr.replace(']', ' ')
        confF.write(writestr)
        confF.write('\nOverall Accuracy: ' + str(allAccuracy))
        confF.write('\n')
        confF.write('Confusion Matrix Top 70%\n')
        writestr = np.array2string(confusionMatrixTop, separator='\t')
        writestr = writestr.replace('[', ' ')
        writestr = writestr.replace(']', ' ')
        confF.write(writestr)
        confF.write('\nAccuracy Top 70%: ' + str(accuracyTop70))
        confF.write('\n')
        confF.write('Confusion Matrix Above 0.90\n')
        writestr = np.array2string(confusionMatrixAbove90, separator='\t')
        writestr = writestr.replace('[', ' ')
        writestr = writestr.replace(']', ' ')
        confF.write(writestr)
        confF.write('\nAccuracy Above 0.90: ' + str(accuracyAbove90))
        confF.write('\n')
        confF.write('Confusion Matrix 0.80_0.90\n')
        writestr = np.array2string(confusionMatrix_b_80_90, separator='\t')
        writestr = writestr.replace('[', ' ')
        writestr = writestr.replace(']', ' ')
        confF.write(writestr)
        confF.write('\nAccuracy 0.80 < Conf <= 0.90: ' + str(accuracy_b_80_90))
        confF.write('\n')
        confF.write('Confusion Matrix 0.70_0.80\n')
        writestr = np.array2string(confusionMatrix_b_70_80, separator='\t')
        writestr = writestr.replace('[', ' ')
        writestr = writestr.replace(']', ' ')
        confF.write(writestr)
        confF.write('\nAccuracy 0.70 < Conf <= 0.80: ' + str(accuracy_b_70_80))
        confF.write('\n')
        confF.write('Confusion Matrix <= 0.70\n')
        writestr = np.array2string(confusionMatrixBelow70, separator='\t')
        writestr = writestr.replace('[', ' ')
        writestr = writestr.replace(']', ' ')
        confF.write(writestr)
        confF.write('\nAccuracy Conf <= 0.70: ' + str(accuracyBelow70))

    allAccuracy, lowerAcc, upperAcc = CM.modelAccuracy(predDF_unsorted['PredictedCluster'], predDF_unsorted['TrueCluster'])
    Kappa, lowerKappa, upperKappa = CM.modelKappa(predDF_unsorted['Confidence'], np.equal(predDF_unsorted['PredictedCluster'], predDF_unsorted['TrueCluster']))
    performance, lowerPerf, upperPerf = CM.modelPerformance(allAccuracy, Kappa, lowerAcc, upperAcc, lowerKappa, upperKappa)

    print(name)
    model_evals[name] = {'Accuracy': allAccuracy, 'Kappa': Kappa, 'Performance': performance,
                         'lowerAcc': lowerAcc, 'upperAcc': upperAcc,
                         'lowerKappa': lowerKappa, 'upperKappa': upperKappa,
                         'lowerPerf': lowerPerf, 'upperPerf': upperPerf}

    predDF.to_csv("../evaluation_test_set/preds_" + name + ext, sep='\t')

evals = pd.DataFrame.from_dict(model_evals, orient='index')
evals.to_csv('../evaluation_test_set/test_set_metric_evaluations.tsv', sep='\t')
