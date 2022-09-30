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
testing_file = '../data_tables/train_test_sets/TestingSet_149Subset_May2021.txt'
qval_file = '../data_tables/qval_dfs/fisher_exact_5x2.Sep_23_2022.combined.tsv'

test_set = list(pd.read_csv(testing_file, sep='\t', header=None, index_col=0).index)
qval_df = pd.read_csv(qval_file, sep='\t', index_col=0)
sig_genes = list(qval_df.loc[qval_df['q'] < 0.10].index)

full_path = '../saved_models/FINALMODEL*/*'
files = glob.glob(full_path)
files.sort(key=lambda x: int(x.split('_')[-1]))
nets = []
for file in files:
    loadedNet = torch.load(file)
    NFEATURES = (len(list(loadedNet.items())[1][1][1]))
    net = nn.Net(10, NFEATURES, 5)
    net.load_state_dict(loadedNet)
    net.eval()
    nets.append(net)


df, targets = format_data.format_inputs(GSM, targetfile, test_set)

df2 = df[df.columns.drop(list(df.filter(regex='.AMP|.DEL|SV')))]
unique, counts = np.unique(df2.values, return_counts=True)
counts = dict(zip(unique, counts))
probabilities_mutations = {1: counts[1] / (counts[1] + counts[2]),
                           2: counts[2]/(counts[1] + counts[2])}

df2 = df[list(df.filter(regex='.AMP|.DEL'))]
unique, counts = np.unique(df2.values, return_counts=True)
counts = dict(zip(unique, counts))
probabilities_scnas = {1: counts[1] / (counts[1] + counts[2]),
                       2: counts[2]/(counts[1] + counts[2])}

results_df = pd.DataFrame(columns=['addInFraction', 'accuracyAll', 'accuracyTop',
                                   'Kappa', 'Performance', 'meanconfidence', 'stdconfidence',
                                   'lowerAccuracy', 'upperAccuracy', 'lowerKappa', 'upperKappa', 'lowerPerformance', 'upperPerformance',
                                   'sv_count_total', 'sv_count_added',
                                   'mut_count_total', 'mut_count_added',
                                   'amp_count_total', 'amp_count_added',
                                   'del_count_total', 'del_count_added'])

df_original = pd.read_csv(GSM, sep='\t', index_col=0)
df_original = df_original.transpose()
df_original = df_original.loc[df_original.index.isin(test_set)]
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

stepsize = 0.05
steps = np.arange(0, 1 + 0.05, stepsize)
pMax = None

for step in steps:
    df = pd.read_csv(GSM, sep='\t', index_col=0)
    df = df.transpose()
    df = df.loc[df.index.isin(test_set)]
    df = df[sig_genes]
    df = df.astype(float)
    df = df.transpose()
    for idx, row in df.iterrows():
        for i in range(0, len(row)):
            curr_col = df.columns[i]
            current_flip = np.random.binomial(1, step, 1)[0]
            if current_flip:
                if 'SV' in idx:
                    val = 3
                elif '.AMP' in idx or '.DEL' in idx:
                    vals_scna = list(probabilities_scnas.keys())
                    probs_scna = list(probabilities_scnas.values())
                    val = np.random.choice(vals_scna, 1, p=probs_scna)[0]
                else:
                    vals_mut = list(probabilities_mutations.keys())
                    probs_mut = list(probabilities_mutations.values())
                    val = np.random.choice(vals_mut, 1, p=probs_mut)[0]
                df.at[idx, curr_col] = val

    mut_count_new = (df.loc[mut_events] != 0).sum().sum()
    mut_count_added = mut_count_new - mut_counts

    amp_count_new = (df.loc[amp_events] != 0).sum().sum()
    amp_count_added = amp_count_new - amp_counts

    del_count_new = (df.loc[del_events] != 0).sum().sum()
    del_count_added = del_count_new - del_counts

    sv_count_new = (df.loc[sv_events] != 0).sum().sum()
    sv_count_added = sv_count_new - sv_counts

    datafile = '../random_add_in_experiment/tmpAddInGSM_fullfeatures_testset.txt'
    df.to_csv(datafile, sep='\t', header=True, index=True)
    data, targets = format_data.format_inputs(datafile, targetfile, test_set,
                                              reduced_version='3.3', drop_empty_vectors=False)

    targets = targets.loc[data.index]
    targets['cluster'] = targets.iloc[:, 0:5].idxmax(axis=1)
    targets['cluster'] = targets['cluster'].map({'C1': 0, 'C2': 1, 'C3': 2, 'C4': 3, 'C5': 4})

    predictionsDF = pd.DataFrame()
    for i in range(len(nets)):
        currNet = nets[i]
        for sample, row in data.iterrows():
            netinputs = torch.tensor(row, dtype=torch.float)
            outputvec = currNet.forward(netinputs).detach().numpy()
            if sample not in predictionsDF.index:
                entry = pd.DataFrame([outputvec], columns=['C1', 'C2', 'C3', 'C4', 'C5'], index=[sample])
                predictionsDF = predictionsDF.append(entry)
            else:
                predictionsDF.loc[sample] += outputvec

    predictionsDF = predictionsDF.div(len(nets))
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

    predictionsDF['Confidence'] = predictionsDF.iloc[:, 0:5].max(axis=1)
    predictionsDF['True Cluster'] = correspondingTargets
    predictionsDF['Predicted Cluster'] = predictionsDF.iloc[:, 0:5].idxmax(axis=1).map({'C1': 0, 'C2': 1, 'C3': 2, 'C4': 3, 'C5': 4})
    predictionsDF = predictionsDF.sort_values(by=['Confidence'])

    meanconfidence = np.mean(predictionsDF['Confidence'])
    stdconfidence = np.std(predictionsDF['Confidence'])

    top70thActuals = predictionsDF.iloc[int(len(predictionsDF.index) * 0.3):]['True Cluster']
    top70thPreds = predictionsDF.iloc[int(len(predictionsDF.index) * 0.3):]['Predicted Cluster']
    top70thAccuracy = np.sum(np.equal(top70thActuals, top70thPreds)) / len(top70thPreds)

    lrx = np.array(predictionsDF['Confidence'])
    lry = np.array(np.equal(predictionsDF['Predicted Cluster'], predictionsDF['True Cluster']).astype(int))

    allAccuracy, lowerAcc, upperAcc = CM.modelAccuracy(predictionsDF['Predicted Cluster'], predictionsDF['True Cluster'])
    Kappa, lowerKappa, upperKappa = CM.modelKappa(lrx, lry)
    performance, lowerPerf, upperPerf = CM.modelPerformance(allAccuracy, Kappa, lowerAcc, upperAcc, lowerKappa, upperKappa)

    print('Frac: ' + str(step), 'Acc:', allAccuracy, 'Kappa:', Kappa, 'Performance: ' + str(performance))
    results_df.loc[step] = [np.round(step, decimals=2), allAccuracy, top70thAccuracy, Kappa, performance,
                            meanconfidence, stdconfidence, lowerAcc, upperAcc,
                            lowerKappa, upperKappa, lowerPerf, upperPerf,
                            sv_count_new, sv_count_added,
                            mut_count_new, mut_count_added,
                            amp_count_new, amp_count_added,
                            del_count_new, del_count_added]

results_df.to_csv('../random_add_in_experiment/resultstable_random_add_in_testset.txt', sep='\t', index=False)
