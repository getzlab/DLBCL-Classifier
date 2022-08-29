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

random.seed(DEFAULT_SEED)
np.random.seed(DEFAULT_SEED)

GSM = "../data_tables/gsm/DLBCL_Staudt_Shipp_CL.for_classifier_training.classifier_subset.fix_sv.fix_ploidy.17-Aug-2022.txt"
targetfile = "../data_tables/confidence_tables/baseline_probabilities.connectivity_based.sensitivity_power2.Aug_17_2022.tsv"
testset_file = '../data_tables/train_test_sets/TestingSet_149Subset_May2021.txt'
samples_file = '../data_tables/sample_sets/ShippStaudtSets.purity0.2.txt'
qval_file = '../data_tables/qval_dfs/fisher_exact_5x2_17-Aug-2022.combined.tsv'
ccf_file = '../data_tables/gsm/DLBCL_700.17-Dec-2021.txt'

ccf_gsm = pd.read_csv(ccf_file, sep='\t', index_col=0)

ccf_gsm.index = ['X' + driver if ('.AMP' in driver or '.DEL' in driver) else driver for driver in ccf_gsm.index]
ccf_gsm.index = ccf_gsm.index.str.replace('_', '.').str.replace('-', '.')

ccf_gsm.loc['MYD88.CCF'] = ccf_gsm.loc[['MYD88.L265P.CCF', 'MYD88.OTHER.CCF']].max(axis=0)

ccf_gsm = ccf_gsm.fillna(0)

# ccf_gsm.loc['SV.BCL2.CCF'] = 1
# ccf_gsm.loc['SV.BCL6.CCF'] = 1
# ccf_gsm.loc['SV.MYC.CCF'] = 1

samples = pd.read_csv(samples_file, sep='\t', index_col=0)
shipp_samples = samples.loc[samples['cohort'] == 'Shipp'].index
testing_set = list(pd.read_csv(testset_file, sep='\t', header=None, index_col=0).index)
testing_set = [x for x in testing_set if x in shipp_samples]
qval_df = pd.read_csv(qval_file, sep='\t', index_col=0)
sig_genes = list(qval_df.loc[qval_df['q'] < 0.10].index)

full_path = '../saved_models/FINALMODEL*/*'
files = glob.glob(full_path)
nets = []

for file in files:
    loadedNet = torch.load(file)
    NFEATURES = (len(list(loadedNet.items())[1][1][1]))
    net = nn.Net(10, NFEATURES, 5)
    net.load_state_dict(loadedNet)
    net.eval()
    nets.append(net)


stepsize = 0.05
steps = np.arange(0, 1 + stepsize, stepsize)

df_original = pd.read_csv(GSM, sep='\t', index_col=0)
df_original = df_original.transpose()
df_original = df_original.loc[df_original.index.isin(testing_set)]
df_original = df_original[sig_genes]
df_original = df_original.astype(float)
df_original = df_original.transpose()

results_df = pd.DataFrame(columns=['ccf_threshold', 'accuracyAll', 'accuracyTop', 'Kappa', 'Performance',
                                   'meanconfidence', 'stdconfidence', 'lowerAccuracy', 'upperAccuracy',
                                   'lowerKappa', 'upperKappa', 'lowerPerformance', 'upperPerformance',
                                   'mut_count_dropped', 'mut_count_kept',
                                   'del_count_dropped', 'del_count_kept',
                                   'amp_count_dropped', 'amp_count_kept',
                                   'sv_count_dropped', 'sv_count_kept'])

ccf_gsm = ccf_gsm.loc[df_original.index + '.CCF', testing_set]

sv_events = df_original.index[df_original.index.str.contains('SV.', regex=False)]
amp_events = df_original.index[df_original.index.str.contains('.AMP', regex=False)]
del_events = df_original.index[df_original.index.str.contains('.DEL', regex=False)]
mut_events = list(set(df_original.index) - set(sv_events) - set(amp_events) - set(del_events))

mut_counts = int((df_original.loc[mut_events] != 0).sum().sum())
del_counts = int((df_original.loc[del_events] != 0).sum().sum())
amp_counts = int((df_original.loc[amp_events] != 0).sum().sum())
sv_counts = int((df_original.loc[sv_events] != 0).sum().sum())

ccf_gsm.loc[[x + '.CCF' for x in mut_events]] = ccf_gsm.loc[[x + '.CCF' for x in mut_events]] + 0.01

pMax = None

for step in steps:
    df = pd.read_csv(GSM, sep='\t', index_col=0)
    df = df.transpose()
    df = df.loc[df.index.isin(testing_set)]
    df = df[sig_genes]
    df = df.astype(float).astype(int)
    df = df.transpose()

    keep_events = (ccf_gsm >= step).astype(int)
    keep_events.index = keep_events.index.str.replace('.CCF', '', regex=False)

    df = df * keep_events

    mut_count_kept = (df.loc[mut_events] != 0).sum().sum()
    mut_count_dropped = mut_counts - mut_count_kept

    amp_count_kept = (df.loc[amp_events] != 0).sum().sum()
    amp_count_dropped = amp_counts - amp_count_kept

    del_count_kept = (df.loc[del_events] != 0).sum().sum()
    del_count_dropped = del_counts - del_count_kept

    sv_count_kept = (df.loc[sv_events] != 0).sum().sum()
    sv_count_dropped = sv_counts - sv_count_kept

    datafile = '../ccf_threshold_experiment/experiment_gsms/thresholded_gsm_thresh' + str(round(step, 2)) + '_test.txt'
    df.to_csv(datafile, sep='\t', header=True, index=True)
    data, targets = format_data.format_inputs(datafile, targetfile, testing_set,
                                              reduced_version='3.2', remove_largest_n=5,
                                              drop_empty_vectors=False)

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
    top70thAccuracy = np.sum(np.equal(top70thActuals, top70thPreds)) / len(top70thActuals)

    lrx = np.array(predictionsDF['Confidence'])
    lry = np.array(np.equal(predictionsDF['Predicted Cluster'], predictionsDF['True Cluster']).astype(int))

    allAccuracy, lowerAcc, upperAcc = CM.modelAccuracy(predictionsDF['Predicted Cluster'], predictionsDF['True Cluster'])
    Kappa, lowerKappa, upperKappa = CM.modelKappa(lrx, lry)
    performance, lowerPerf, upperPerf = CM.modelPerformance(allAccuracy, Kappa, lowerAcc, upperAcc, lowerKappa, upperKappa)

    print('Threshold: ' + str(np.round(step, decimals=2)), 'Acc:', allAccuracy, 'Kappa:', Kappa, 'Performance: ' + str(performance))
    results_df.loc[step] = [str(np.round(step, decimals=2)), allAccuracy, top70thAccuracy, Kappa, performance,
                            meanconfidence, stdconfidence, lowerAcc, upperAcc,
                            lowerKappa, upperKappa, lowerPerf, upperPerf,
                            mut_count_dropped, mut_count_kept,
                            del_count_dropped, del_count_kept,
                            amp_count_dropped, amp_count_kept,
                            sv_count_dropped, sv_count_kept]

results_df.to_csv('../ccf_threshold_experiment/resultstable_ccf_threshold_test.txt', sep='\t', index=False)
