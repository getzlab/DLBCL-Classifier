import glob
import pandas as pd
import numpy as np
import torch
import nn
import random
import format_data
import calculate_model_metrics as CM
from matplotlib import pyplot as plt

DEFAULT_SEED = 123
WINDOW_SIZE = 10
STEP_SIZE = 10
NUM_CORRELATION_ITERS = 1000

random.seed(DEFAULT_SEED)
np.random.seed(DEFAULT_SEED)

CLUSTER_LABELS = ['C' + str(x) for x in range(1, 6)]

GSM = "../data_tables/gsm/DLBCL.699.fullGSM.Sep_23_2022.tsv"
targetfile = "../data_tables/confidence_tables/baseline_probabilities.connectivity_based.sensitivity_power2.Sep_23_2022.tsv"
training_file = '../data_tables/train_test_sets/TrainingSet_550Subset_May2021.txt'
qval_file = '../data_tables/qval_dfs/fisher_exact_5x2.Sep_23_2022.combined.tsv'
purity_file = '../data_tables/purities_ploidies/ALLPurities_fixednames.tsv'
ploidy_file = '../data_tables/purities_ploidies/PloidyDataFrame.txt'
sample_set_file = '../data_tables/sample_sets/ShippStaudtSets.purity0.2.txt'
preds_file = '../evaluation_validation_set/confidence_adjusted_tables/NN_reducedV3.3_nfeatures21_pMax0.93344957.tsv'
alt_counts_file = '../reduce_purity_experiment/gsm_alt_counts.tsv'
sv_alt_counts_file = '../reduce_purity_experiment/sv_alt_counts.tsv'
cn_file = '../reduce_purity_experiment/copy_ratio_mat.tsv'


alt_counts = pd.read_csv(alt_counts_file, sep='\t', index_col=0)
alt_counts.index = alt_counts.index.str.replace('-', '.')

sv_alt_counts = pd.read_csv(sv_alt_counts_file, sep='\t', index_col=0)

alt_counts = pd.concat([alt_counts, sv_alt_counts])

cr_0s = pd.read_csv(cn_file, sep='\t', index_col=0, low_memory=False)
cr_0s.index = 'X' + cr_0s.index.str.replace(':', '.').str.replace('_', '.').str.upper()
cr_0s = cr_0s + 1

purities = pd.read_csv(purity_file, sep='\t', index_col=0, header=None)
purities.columns = ['purity']
ploidies = pd.read_csv(ploidy_file, sep='\t', index_col=0, header=None)
ploidies.columns = ['ploidy']

purities = purities.astype(float)
ploidies = ploidies.astype(float)

labels = pd.read_csv(targetfile, sep='\t', index_col=0)
labels['cluster'] = labels['cluster'].map({1: 'C1', 2: 'C2', 3: 'C3', 4: 'C4', 5: 'C5'})

sample_sets = pd.read_csv(sample_set_file, index_col=0, sep='\t')
shipp_samples = sample_sets.loc[sample_sets['cohort'] == 'Shipp'].index

# DLBCL_C_D_1148_NULLPAIR missing tumor read counts for SV event. exclude sample
shipp_samples = [x for x in shipp_samples if x != 'DLBCL_C_D_1148_NULLPAIR']

preds = pd.read_csv(preds_file, sep='\t', index_col=0)
preds.columns = CLUSTER_LABELS
preds['cluster'] = preds.idxmax(axis=1)

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
    validationfile = '../all_validation_sets/NN_evaluation_seeds1_100_folds5_reducedV3.3/NN_evaluation_seeds1_100_folds5_reducedV3.3' \
                     + '_' + str(curriter + 1) + '_' + str(currfold)
    validationSamples = list(pd.read_csv(validationfile, sep='\t', index_col=0, header=None).index)
    validation_sets[netnum] = validationSamples
    loadedNet = torch.load(file)
    NFEATURES = (len(list(loadedNet.items())[1][1][1]))
    net = nn.Net(10, NFEATURES, 5)
    net.load_state_dict(loadedNet)
    net.eval()
    nets[netnum] = net

# starting_samples = {x: None for x in CLUSTER_LABELS}
# all_starting_samples = []
preds['purity'] = purities['purity']
labels['purity'] = purities['purity']

all_starting_samples = purities.loc[purities.index.isin(shipp_samples) & (purities['purity'] >= 0.50) & (purities.index.isin(preds.index))].index

preds = preds.sort_values(by=['cluster', 'purity'], ascending=[True, False])
labels = labels.sort_values(by=['cluster', 'purity'], ascending=[True, False])

# n_per_clus = 10
# for clus in CLUSTER_LABELS:
#     labels_c = labels.loc[(labels['cluster'] == clus) & (labels.index.isin(shipp_samples)) & (labels.index.isin(preds.index))]
#     starting_samples[clus] = list(labels_c.iloc[0:n_per_clus].index)
#     all_starting_samples = all_starting_samples + list(labels_c.iloc[0:n_per_clus].index)

df_original = pd.read_csv(GSM, sep='\t', index_col=0)
df_original = df_original.transpose()
df_original = df_original.loc[df_original.index.isin(all_starting_samples)]
df_original = df_original[sig_genes]
df_original = df_original.astype(float)
df_original = df_original.transpose()

results_df = pd.DataFrame(index=df_original.columns)
alt_counts = alt_counts.loc[alt_counts.index.isin(df_original.index)]

pMax = None


def power_calc(CR_0s, purity_sim):
    event = CR_0s.name
    samples = CR_0s.index
    purities_vec = purities.loc[samples, 'purity']
    ploidies_vec = ploidies.loc[samples, 'ploidy']

    Q_vec = (CR_0s * (purities_vec * ploidies_vec + 2*(1-purities_vec)) - 2*(1-purities_vec)) / purities_vec

    CR_vec = (purity_sim*Q_vec+2*(1-purity_sim)) / (purity_sim*ploidies_vec+2*(1-purity_sim))

    new_vals = CR_vec.copy(deep=True)
    if '.AMP' in event:
        for n, e in CR_vec.items():
            if e >= 1.9:
                new_vals[n] = 2
            elif e < 1.9 and e >= 1.1:
                new_vals[n] = 1
            else:
                new_vals[n] = 0
    elif '.DEL' in event:
        for n, e in CR_vec.items():
            if e < 0.1:
                new_vals[n] = 2
            elif e < 0.9 and e >= 0.1:
                new_vals[n] = 1
            else:
                new_vals[n] = 0
    else:
        print('this should never happen')
        exit()

    return new_vals


steps = [0.50, 0.40, 0.30, 0.20, 0.15, 0.10, 0.05, 0.02]
cnas = [x for x in df_original.index if '.AMP' in x or '.DEL' in x]
cr_0s = cr_0s[all_starting_samples]

results_df['initial_pred'] = preds['cluster'].map({'C1': 1, 'C2': 2, 'C3': 3, 'C4': 4, 'C5': 5})
results_df['initial_confidence'] = preds.iloc[:, 0:5].max(axis=1)

cr_0s = cr_0s.loc[cr_0s.index.isin(sig_genes)]

amp_cnas = [x for x in df_original.index if '.AMP' in x]
del_cnas = [x for x in df_original.index if '.DEL' in x]
svs = [x for x in df_original.index if 'SV.' in x]
muts = [x for x in df_original.index if x not in amp_cnas and x not in del_cnas and x not in svs]

count_df = pd.DataFrame(columns=['initial_count'] + ['count_' + str(x) for x in steps], index=['mut_nonsil', 'mut_sil',
                                                                                               'sv',
                                                                                               'amp_2', 'amp_1',
                                                                                               'del_2', 'del_1'])

count_df.loc['mut_nonsil', 'initial_count'] = (df_original.loc[muts].astype(int) == 2).sum(axis=1).sum()
count_df.loc['mut_sil', 'initial_count'] = (df_original.loc[muts].astype(int) == 1).sum(axis=1).sum()
count_df.loc['sv', 'initial_count'] = (df_original.loc[svs].astype(int) == 3).sum(axis=1).sum()
count_df.loc['amp_2', 'initial_count'] = (df_original.loc[amp_cnas].astype(int) == 2).sum(axis=1).sum()
count_df.loc['amp_1', 'initial_count'] = (df_original.loc[amp_cnas].astype(int) == 1).sum(axis=1).sum()
count_df.loc['del_2', 'initial_count'] = (df_original.loc[del_cnas].astype(int) == 2).sum(axis=1).sum()
count_df.loc['del_1', 'initial_count'] = (df_original.loc[del_cnas].astype(int) == 1).sum(axis=1).sum()

for step in steps:
    curr_purities = purities.loc[purities.index.isin(all_starting_samples)].copy(deep=True)
    cnas_orig = df_original.loc[df_original.index.str.contains('.AMP') | df_original.index.str.contains('.DEL'), all_starting_samples]

    cr_step_mask = cr_0s.apply(lambda x: power_calc(x, step), axis=1)

    alt_counts_step = alt_counts[all_starting_samples].copy(deep=True).apply(lambda x: np.floor(x * (step / curr_purities.loc[x.name, 'purity'])))
    alt_counts_mask = alt_counts_step > 2

    df_power = df_original.copy(deep=True).loc[alt_counts_mask.index, all_starting_samples]
    df_power = df_power * alt_counts_mask

    df_power = pd.concat([df_power, cr_step_mask], axis=0)

    count_df.loc['mut_nonsil', 'count_' + str(step)] = (df_power.loc[muts].astype(int) == 2).sum(axis=1).sum()
    count_df.loc['mut_sil', 'count_' + str(step)] = (df_power.loc[muts].astype(int) == 1).sum(axis=1).sum()
    count_df.loc['sv', 'count_' + str(step)] = (df_power.loc[svs].astype(int) == 3).sum(axis=1).sum()
    count_df.loc['amp_2', 'count_' + str(step)] = (df_power.loc[amp_cnas].astype(int) == 2).sum(axis=1).sum()
    count_df.loc['amp_1', 'count_' + str(step)] = (df_power.loc[amp_cnas].astype(int) == 1).sum(axis=1).sum()
    count_df.loc['del_2', 'count_' + str(step)] = (df_power.loc[del_cnas].astype(int) == 2).sum(axis=1).sum()
    count_df.loc['del_1', 'count_' + str(step)] = (df_power.loc[del_cnas].astype(int) == 1).sum(axis=1).sum()

    datafile = '../reduce_purity_experiment/experiment_gsms/reducepurity_fullfeatures_step' + str(step) + '.txt'
    df_power.to_csv(datafile, sep='\t', header=True, index=True)
    data, targets = format_data.format_inputs(datafile, targetfile, all_starting_samples,
                                              reduced_version='3.3',
                                              drop_empty_vectors=False)

    targets = targets.loc[data.index]
    targets['cluster'] = targets.iloc[:, 0:5].idxmax(axis=1)
    targets['cluster'] = targets['cluster'].map({'C1': 0, 'C2': 1, 'C3': 2, 'C4': 3, 'C5': 4})

    predictionsDF = pd.DataFrame()
    for i in range(len(nets)):
        currNet = nets[i]
        currVal = validation_sets[i]
        valDF = data.loc[[x for x in currVal if x in all_starting_samples]]
        for sample, row in valDF.iterrows():
            netinputs = torch.tensor(row, dtype=torch.float)
            outputvec = currNet.forward(netinputs).detach().numpy()
            if sample not in predictionsDF.index:
                entry = pd.DataFrame([outputvec], columns=['C1', 'C2', 'C3', 'C4', 'C5'], index=[sample])
                predictionsDF = predictionsDF.append(entry)
            else:
                predictionsDF.loc[sample] += outputvec

    predictionsDF = predictionsDF.div(100)

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

    predictionsDF.to_csv('../reduce_purity_experiment/preds/reducepurity_preds_' + '_step' + str(step) + '.txt', sep='\t')

    results_df['clus' + str(step)] = predictionsDF['Predicted Cluster'] + 1
    results_df['conf' + str(step)] = predictionsDF['Confidence']

results_df.to_csv('../reduce_purity_experiment/resultstable_reducepurity.txt', sep='\t')
count_df.to_csv('../reduce_purity_Experiment/countstable_reducepurity.txt', sep='\t')
