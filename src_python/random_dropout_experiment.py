import glob
import pandas as pd
import numpy as np
import torch
import nn
import random
import format_data
import calculate_model_metrics as CM

DEFAULT_SEED = 123
WINDOW_SIZE = 10
STEP_SIZE = 10
NUM_CORRELATION_ITERS = 1000
TEST_MODE = False

random.seed(DEFAULT_SEED)
np.random.seed(DEFAULT_SEED)

GSM = '../data_tables/gsm/DLBCL.699.fullGSM.Sep_23_2022.tsv'
targetfile = '../data_tables/confidence_tables/baseline_probabilities.connectivity_based.sensitivity_power2.Sep_23_2022.tsv'
training_file = '../data_tables/train_test_sets/TrainingSet_550Subset_May2021.txt'
qval_file = '../data_tables/qval_dfs/fisher_exact_5x2.Sep_23_2022.combined.tsv'

training_set = list(pd.read_csv(training_file, sep='\t', header=None, index_col=0).index)
qval_df = pd.read_csv(qval_file, sep='\t', index_col=0)
sig_genes = list(qval_df.loc[qval_df['q'] < 0.10].index)

full_path = '../saved_models/FINALMODEL*/*'
files = glob.glob(full_path)
files.sort(key=lambda x: int(x.split('_')[-1]))
nets = [None] * len(files)
validation_sets = [None] * len(files)
for file in files:
    netnum = int(file.split('_')[-1]) - 1
    curriter = int(np.floor(netnum/5))
    currfold = (netnum % 5)
    validationfile = '../all_validation_sets/NN_evaluation_seeds1_100_folds5_reducedV3.4_removeN5/NN_evaluation_seeds1_100_folds5_reducedV3.4_removeN5' \
                     + '_' + str(curriter + 1) + '_' + str(currfold)
    validationSamples = list(pd.read_csv(validationfile, sep='\t', index_col=0, header=None).index)
    validation_sets[netnum] = validationSamples
    loadedNet = torch.load(file)
    NFEATURES = (len(list(loadedNet.items())[1][1][1]))
    net = nn.Net(10, NFEATURES, 5)
    net.load_state_dict(loadedNet)
    net.eval()
    nets[netnum] = net

stepsize = 0.05
steps = np.arange(0.05, 1 + stepsize, stepsize)
steps = steps[::-1]

pMax = None

if TEST_MODE:
    fn = '../random_dropout_experiment/resultstable_random_dropout_test.txt'
else:
    fn = '../random_dropout_experiment/resultstable_random_dropout.txt'


df_original = pd.read_csv(GSM, sep='\t', index_col=0)
df_original = df_original.transpose()
df_original = df_original.loc[df_original.index.isin(training_set)]
df_original = df_original[sig_genes]
df_original = df_original.astype(float)
df_original = df_original.transpose()

sv_events = df_original.index[df_original.index.str.contains('SV.', regex=False)]
amp_events = df_original.index[df_original.index.str.contains('.AMP', regex=False)]
del_events = df_original.index[df_original.index.str.contains('.DEL', regex=False)]
mut_events = list(set(df_original.index) - set(sv_events) - set(amp_events) - set(del_events))

mut_counts = int((df_original.loc[mut_events] != 0).sum().sum())
del_counts = int((df_original.loc[del_events] != 0).sum().sum())
amp_counts = int((df_original.loc[amp_events] != 0).sum().sum())
sv_counts = int((df_original.loc[sv_events] != 0).sum().sum())

df_original = df_original.transpose()

df_original = format_data.construct_reduced_winning_version(df_original)

masks = []
for step in steps:
    mask = np.random.binomial(1, step, (163, 550))
    masks.append(mask)

results_df = pd.DataFrame(columns=['dropout_probability', 'accuracyAll', 'accuracyTop',
                                   'Kappa', 'Performance', 'meanconfidence', 'stdconfidence',
                                   'lowerAccuracy', 'upperAccuracy', 'lowerKappa', 'upperKappa', 'lowerPerformance', 'upperPerformance',
                                   'sv_count_total', 'sv_count_dropped',
                                   'mut_count_total', 'mut_count_dropped',
                                   'amp_count_total', 'amp_count_dropped',
                                   'del_count_total', 'del_count_dropped'])

for step in steps:
    df = pd.read_csv(GSM, sep='\t', index_col=0)
    df = df.transpose()
    df = df.loc[df.index.isin(training_set)]
    df = df[sig_genes]
    df = df.astype(float).astype(int)
    df = df.transpose()
    mask = masks.pop(0)
    df = df * mask
    # if step < 0.10:
    #     random_events = np.random.binomial(1, 0.01, (df.shape[0], df.shape[1]))
    #     df = df + random_events
    datafile = '../random_dropout_experiment/experiment_gsms/tmpRandomDroppedGSM_fullfeatures_step' + str(round(1-step, 2)) + '.txt'
    df.to_csv(datafile, sep='\t', header=True, index=True)
    data, targets = format_data.format_inputs(datafile, targetfile, training_set,
                                              reduced_version='3.4',
                                              remove_largest_n=5,
                                              drop_empty_vectors=False)

    mut_count_new = (df.loc[mut_events] != 0).sum().sum()
    mut_count_dropped = mut_counts - mut_count_new

    amp_count_new = (df.loc[amp_events] != 0).sum().sum()
    amp_count_dropped = amp_counts - amp_count_new

    del_count_new = (df.loc[del_events] != 0).sum().sum()
    del_count_dropped = del_counts - del_count_new

    sv_count_new = (df.loc[sv_events] != 0).sum().sum()
    sv_count_dropped = sv_counts - sv_count_new

    targets = targets.loc[data.index]
    targets['cluster'] = targets.iloc[:, 0:5].idxmax(axis=1)
    targets['cluster'] = targets['cluster'].map({'C1': 0, 'C2': 1, 'C3': 2, 'C4': 3, 'C5': 4})

    predictionsDF = pd.DataFrame()
    for i in range(len(nets)):
        currNet = nets[i]
        currVal = validation_sets[i]
        valDF = data.loc[currVal]
        for sample, row in valDF.iterrows():
            netinputs = torch.tensor(row, dtype=torch.float)
            outputvec = currNet.forward(netinputs).detach().numpy()
            if sample not in predictionsDF.index:
                entry = pd.DataFrame([outputvec], columns=['C1', 'C2', 'C3', 'C4', 'C5'], index=[sample])
                predictionsDF = predictionsDF.append(entry)
            else:
                predictionsDF.loc[sample] += outputvec

    predictionsDF = predictionsDF.div(100)
    if TEST_MODE:
        predictionsDF.to_csv('../random_dropout_experiment/predictions_raw_test_step' + str(1-step) + '.txt')
    if not pMax:
        pMax = predictionsDF.max().max()
        print('Setting pMax', pMax)

    correspondingTargets = []
    for idx, row in predictionsDF.iterrows():
        newArr = np.float64(np.array(row))
        pWin = np.float64(np.max(row))
        newConf = np.float64(pWin / pMax)
        shrinkfactor = np.float64((pMax - pWin) / (pMax * np.float64(np.sum(np.delete(newArr, newArr.argmax())))))
        for idx2 in range(len(newArr)):
            if idx2 == newArr.argmax():
                newArr[idx2] = newConf
            else:
                newArr[idx2] = newArr[idx2] * shrinkfactor
        predictionsDF.loc[idx] = newArr
        correspondingTargets.append(targets.loc[idx]['cluster'])

    if TEST_MODE:
        predictionsDF.to_csv('../random_dropout_experiment/predictions_scaled_test_step' + str(1-step) + '.txt')

    predictionsDF['Confidence'] = predictionsDF.iloc[:, 0:5].max(axis=1)
    predictionsDF['True Cluster'] = correspondingTargets
    predictionsDF['Predicted Cluster'] = predictionsDF.iloc[:, 0:5].idxmax(axis=1).map({'C1': 0, 'C2': 1, 'C3': 2, 'C4': 3, 'C5': 4})
    predictionsDF = predictionsDF.sort_values(by=['Confidence'])

    meanconfidence = np.mean(predictionsDF['Confidence'])
    stdconfidence = np.std(predictionsDF['Confidence'])

    top70thActuals = predictionsDF.iloc[int(len(predictionsDF.index) * 0.3):]['True Cluster']
    top70thPreds = predictionsDF.iloc[int(len(predictionsDF.index) * 0.3):]['Predicted Cluster']
    top70thAccuracy = np.sum(np.equal(top70thActuals, top70thPreds)) / len(top70thActuals)

    lrx = np.array(predictionsDF['Confidence'])
    lry = np.array(np.equal(predictionsDF['Predicted Cluster'], predictionsDF['True Cluster']).astype(int))

    allAccuracy, lowerAcc, upperAcc = CM.modelAccuracy(predictionsDF['Predicted Cluster'], predictionsDF['True Cluster'])
    Kappa, lowerKappa, upperKappa = CM.modelKappa(lrx, lry)
    performance, lowerPerf, upperPerf = CM.modelPerformance(allAccuracy, Kappa, lowerAcc, upperAcc, lowerKappa, upperKappa)

    predictionsDF.to_csv('../random_dropout_experiment/preds/dropout_preds_' + str(round(1-step, 2)) + '.txt', sep='\t')

    print('Frac: ' + str(np.round(1 - step, decimals=4)), 'Acc:', allAccuracy, 'Kappa:', Kappa, 'Performance: ' + str(performance))
    results_df.loc[step] = [np.round(1 - step, decimals=2), allAccuracy, top70thAccuracy, Kappa, performance,
                            meanconfidence, stdconfidence, lowerAcc, upperAcc,
                            lowerKappa, upperKappa, lowerPerf, upperPerf,
                            sv_count_new, sv_count_dropped,
                            mut_count_new, mut_count_dropped,
                            amp_count_new, amp_count_dropped,
                            del_count_new, del_count_dropped]

results_df.to_csv('../random_dropout_experiment/resultstable_random_dropout.txt', sep='\t', index=False)
