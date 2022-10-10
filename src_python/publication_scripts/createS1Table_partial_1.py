import pandas as pd
import numpy as np

natmed_tableS1 = pd.read_csv('../../data_tables/2017_11_14_Table_S1.txt',
                             sep='\t', index_col=0, header=0, skiprows=1)
gsm = pd.read_csv("../../data_tables/gsm/DLBCL.699.fullGSM.Sep_23_2022.tsv",
                  sep='\t', index_col=0)
qval_df = pd.read_csv('../../data_tables/qval_dfs/fisher_exact_5x2.Sep_23_2022.combined.tsv',
                      sep='\t', low_memory=False, index_col=0)
coverages = pd.read_csv('../../data_tables/clustering_labels/per_sample_tumor_coverage_formatted.txt',
                        sep='\t', index_col=0)
cootable = pd.read_csv('../../data_tables/survival_data/COOtable.txt',
                       sep='\t', index_col=0)
purities = pd.read_csv('../../data_tables/purities_ploidies/ALLPurities_fixednames.tsv',
                       sep='\t', index_col=0, header=None)
ploidies = pd.read_csv('../../data_tables/purities_ploidies/PloidyDataFrame.txt',
                       sep='\t', index_col=0, header=None)
ploidies.columns = ['PLOIDY']
labels = pd.read_csv('../../data_tables/confidence_tables/baseline_probabilities.connectivity_based.sensitivity_power2.Sep_23_2022.tsv',
                     sep='\t', index_col=0)
nci_samples = list(pd.read_csv('../../data_tables/purities_ploidies/NCIsamples.txt',
                               sep='\t', index_col=0).index)
tcga_samples = list(pd.read_csv('../../data_tables/purities_ploidies/TCGAsamples.txt',
                                sep='\t', index_col=0).index.str.replace('-', '_'))
training_set = list(pd.read_csv('../../data_tables/train_test_sets/TrainingSet_550Subset_May2021.txt',
                                sep='\t', index_col=0, header=None).index)
test_set = list(pd.read_csv('../../data_tables/train_test_sets/TestingSet_149Subset_May2021.txt',
                            sep='\t', index_col=0, header=None).index)

train_preds = pd.read_csv('../../evaluation_validation_set/confidence_adjusted_tables/NN_reducedV3.3_nfeatures21_pMax0.93344957.tsv',
                          sep='\t', index_col=0)

test_preds = pd.read_csv('../../evaluation_test_set/NN_reducedV3.3_nfeatures21_testsetEval.tsv',
                          sep='\t', index_col=0)

raw_counts = pd.read_csv('../../data_tables/raw_muts_cnas_counts.tsv', sep='\t', index_col=0)


staudt_only_labels = pd.read_csv('../../data_tables/clustering_labels/GSM699_cluster_Sep_23_2022_STAUDT.bestclus.txt', sep='\t', index_col=0)
shipp_only_labels = pd.read_csv('../../data_tables/clustering_labels/GSM699_cluster_Sep_23_2022_SHIPP.bestclus.remapped.txt', sep='\t', index_col=0)

train_preds['Confidence'] = train_preds.max(axis=1)
train_preds['PredictedCluster'] = train_preds.iloc[:, 0:5].idxmax(axis=1).astype(int) + 1

train_preds.columns = ['Predicted C1', 'Predicted C2', 'Predicted C3', 'Predicted C4', 'Predicted C5'] + list(train_preds.columns[5::])
test_preds.columns = ['Predicted C1', 'Predicted C2', 'Predicted C3', 'Predicted C4', 'Predicted C5'] + list(test_preds.columns[5::])

test_preds['Above90'] = test_preds['Confidence'] > 0.90
train_preds['Above90'] = train_preds['Confidence'] > 0.90

test_preds['Top70 Perc. Confident'] = test_preds['Confidence'] > np.percentile(test_preds['Confidence'], 30)
train_preds['Top70 Perc. Confident'] = train_preds['Confidence'] > np.percentile(train_preds['Confidence'], 30)


natmed_tableS1.index = natmed_tableS1.index.str.replace('-', '_').str.upper()
natmed_tableS1['pair_id'] = natmed_tableS1['pair_id'].str.replace('-', '_').str.upper().str.replace('_TP.+', '', regex=True)


purities.index = purities.index.str.replace('-', '_')
purities.columns = ['Purity']

labels['confidence'] = labels.max(axis=1)
labels['cluster'] = labels.iloc[:, 0:5].idxmax(axis=1).map({'C1': 1, 'C2': 2, 'C3': 3, 'C4': 4, 'C5': 5})

gsm = gsm.T
sampleset = set(gsm.index) - set(coverages.index)
gsm = gsm.drop(sampleset, axis=0)
todrop = set(gsm.columns) - (set(qval_df.loc[qval_df['q'] < 0.10].index) - {'OR51B6', 'OR10V1'})
gsm = gsm.drop(list(todrop), axis=1)

columns = ['Cohort', 'COO', 'Gender', 'Age-at first diagnosis',
           'Sample Preparation', 'Mean Target Coverage', 'Median Target Coverage',
           'Ploidy (ABSOLUTE)', 'Purity (ABSOLUTE)', 'Pair Status (Tumor Only/Normal)',
           'Number of Mutations', 'Fraction Genome Deleted', 'Fraction Genome Amplified',
           'Number of Drivers - Mutations', 'Number of Drivers - Non-Silent Mutations',
           'Number of Drivers - SCNAs', 'Number of Drivers - Amplifications', 'Number of Drivers - Deletions',
           'Number of Drivers - SVs',
           'Number of WT Mutations (0)', 'Number of WT SCNAs (0)', 'Number of WT Deletions (0)', 'Number of WT Amplifications (0)',
           'Number of WT SVs (0)',
           'Test/Train Set Membership', 'Cluster',
           'Target C1', 'Target C2', 'Target C3', 'Target C4', 'Target C5',
           'Predicted C1', 'Predicted C2', 'Predicted C3', 'Predicted C4', 'Predicted C5', 'PredictedCluster', 'Confidence',
           'Above90', 'Top70 Perc. Confident',
           'Staudt Only Cluster', 'Shipp Only Cluster']

mastertable = pd.DataFrame(np.nan, index=coverages.index, columns=columns)
mastertable.index.name = 'ID'


for idx,row in mastertable.iterrows():
    # Cohort
    if idx in nci_samples:
        mastertable.loc[idx, 'Cohort'] = 'Schmitz et al.'
    elif idx in tcga_samples:
        mastertable.loc[idx, 'Cohort'] = 'Schmitz et al.'
    else:
        mastertable.loc[idx, 'Cohort'] = 'Chapuy et al.'

    # COO
    if idx in cootable.index:
        if cootable.loc[idx, 'COO'] != 'na':
            mastertable.loc[idx, 'COO'] = cootable.loc[idx, 'COO']

    # Misc stats
    if idx in list(natmed_tableS1['pair_id']):

        tmp = natmed_tableS1.loc[natmed_tableS1['pair_id'] == idx, 'Gender'].values[0]
        if tmp != 'na':
            mastertable.loc[idx, 'Gender'] = tmp

        tmp = natmed_tableS1.loc[natmed_tableS1['pair_id'] == idx, 'Age-at first diagnosis'].values[0]
        if tmp != 'na':
            mastertable.loc[idx, 'Age-at first diagnosis'] = tmp

    # Sample prep
    if idx in natmed_tableS1.index:
        mastertable.loc[idx, 'Sample Preparation'] = natmed_tableS1.loc[idx, 'sample preparation']
    elif idx in list(natmed_tableS1['pair_id']):
        tmp = natmed_tableS1.loc[natmed_tableS1['pair_id'] == idx, 'sample preparation'].values[0]
        mastertable.loc[idx, 'Sample Preparation'] = tmp
    elif idx in nci_samples or idx in tcga_samples:
        mastertable.loc[idx, 'Sample Preparation'] = 'Frozen'

    # Ploidy
    if idx in ploidies.index:
        mastertable.loc[idx, 'Ploidy (ABSOLUTE)'] = ploidies.loc[idx, 'PLOIDY']

    # Purity
    if idx in purities.index:
        mastertable.loc[idx, 'Purity (ABSOLUTE)'] = purities.loc[idx, 'Purity']

    # TN/TO
    if idx in nci_samples:
        mastertable.loc[idx, 'Pair Status (Tumor Only/Normal)'] = 'TO'
    elif idx in tcga_samples:
        mastertable.loc[idx, 'Pair Status (Tumor Only/Normal)'] = 'TN'
    elif idx in natmed_tableS1.index:
        mastertable.loc[idx, 'Pair Status (Tumor Only/Normal)'] = natmed_tableS1.loc[idx, 'pair_tumor_only_status']
    elif idx in list(natmed_tableS1['pair_id']):
        tmp = natmed_tableS1.loc[natmed_tableS1['pair_id'] == idx, 'pair_tumor_only_status'].values[0]
        mastertable.loc[idx, 'Pair Status (Tumor Only/Normal)'] = tmp

    if idx in gsm.index:
        muts = ~gsm.columns.str.contains('\.AMP|\.DEL|SV\.')
        cnas = gsm.columns.str.contains('\.AMP|\.DEL')
        dels = gsm.columns.str.contains('\.DEL')
        amps = gsm.columns.str.contains('\.AMP')
        svs = gsm.columns.str.contains('SV\.')

        mutGSM = gsm.loc[:, muts].astype(float).astype(int)
        cnaGSM = gsm.loc[:, cnas].astype(float).astype(int)
        svGSM = gsm.loc[:, svs].astype(float).astype(int)
        delGSM = gsm.loc[:, dels].astype(float).astype(int)
        ampGSM = gsm.loc[:, amps].astype(float).astype(int)

        nummut = (mutGSM.loc[idx] != 0).sum()
        numsilent = (mutGSM.loc[idx] == 1).sum()
        numnonsilent = nummut - numsilent
        numcna = (cnaGSM.loc[idx] != 0).sum()
        numsv = (svGSM.loc[idx] != 0).sum()
        numdel = (delGSM.loc[idx] != 0).sum()
        numamp = (ampGSM.loc[idx] != 0).sum()

        num_m_wt = (mutGSM.loc[idx].astype(float) == 0.0).sum()
        num_sc_wt = (cnaGSM.loc[idx].astype(float) == 0.0).sum()
        num_sv_wt = (svGSM.loc[idx].astype(float) == 0.0).sum()
        num_amp_wt = (ampGSM.loc[idx].astype(float) == 0.0).sum()
        num_del_wt = (delGSM.loc[idx].astype(float) == 0.0).sum()

        mastertable.loc[idx, 'Number of Drivers - Mutations'] = nummut
        mastertable.loc[idx, 'Number of Drivers - Non-Silent Mutations'] = numnonsilent
        mastertable.loc[idx, 'Number of Drivers - SCNAs'] = numcna
        mastertable.loc[idx, 'Number of Drivers - SVs'] = numsv
        mastertable.loc[idx, 'Number of Drivers - Amplifications'] = numamp
        mastertable.loc[idx, 'Number of Drivers - Deletions'] = numdel
        mastertable.loc[idx, 'Number of WT Mutations (0)'] = num_m_wt
        mastertable.loc[idx, 'Number of WT SCNAs (0)'] = num_sc_wt
        mastertable.loc[idx, 'Number of WT Deletions (0)'] = num_del_wt
        mastertable.loc[idx, 'Number of WT Amplifications (0)'] = num_amp_wt
        mastertable.loc[idx, 'Number of WT SVs (0)'] = num_sv_wt

        # 'Number of Mutations', 'Fraction Genome Deleted', 'Fraction Genome Amplified',

        mastertable.loc[idx, 'Number of Mutations'] = raw_counts.loc[idx, 'muts']
        mastertable.loc[idx, 'Fraction Genome Deleted'] = raw_counts.loc[idx, 'del_frac']
        mastertable.loc[idx, 'Fraction Genome Amplified'] = raw_counts.loc[idx, 'amp_frac']

    if idx in training_set and idx in gsm.index:
        mastertable.loc[idx, 'Test/Train Set Membership'] = 'Train'
    elif idx in test_set and idx in gsm.index:
        mastertable.loc[idx, 'Test/Train Set Membership'] = 'Test'

    if idx in labels.index:
        mastertable.loc[idx, 'Cluster'] = labels.loc[idx, 'cluster']
        mastertable.loc[idx, 'Target C1'] = labels.loc[idx, 'C1']
        mastertable.loc[idx, 'Target C2'] = labels.loc[idx, 'C2']
        mastertable.loc[idx, 'Target C3'] = labels.loc[idx, 'C3']
        mastertable.loc[idx, 'Target C4'] = labels.loc[idx, 'C4']
        mastertable.loc[idx, 'Target C5'] = labels.loc[idx, 'C5']

    if idx in training_set:
        mastertable.loc[idx, 'Predicted C1'] = train_preds.loc[idx, 'Predicted C1']
        mastertable.loc[idx, 'Predicted C2'] = train_preds.loc[idx, 'Predicted C2']
        mastertable.loc[idx, 'Predicted C3'] = train_preds.loc[idx, 'Predicted C3']
        mastertable.loc[idx, 'Predicted C4'] = train_preds.loc[idx, 'Predicted C4']
        mastertable.loc[idx, 'Predicted C5'] = train_preds.loc[idx, 'Predicted C5']
        mastertable.loc[idx, 'Confidence'] = train_preds.loc[idx, 'Confidence']
        mastertable.loc[idx, 'PredictedCluster'] = train_preds.loc[idx, 'PredictedCluster']
        mastertable.loc[idx, 'Above90'] = train_preds.loc[idx, 'Above90']
        mastertable.loc[idx, 'Top70 Perc. Confident'] = train_preds.loc[idx, 'Top70 Perc. Confident']

    if idx in test_set:
        mastertable.loc[idx, 'Predicted C1'] = test_preds.loc[idx, 'Predicted C1']
        mastertable.loc[idx, 'Predicted C2'] = test_preds.loc[idx, 'Predicted C2']
        mastertable.loc[idx, 'Predicted C3'] = test_preds.loc[idx, 'Predicted C3']
        mastertable.loc[idx, 'Predicted C4'] = test_preds.loc[idx, 'Predicted C4']
        mastertable.loc[idx, 'Predicted C5'] = test_preds.loc[idx, 'Predicted C5']
        mastertable.loc[idx, 'Confidence'] = test_preds.loc[idx, 'Confidence']
        mastertable.loc[idx, 'PredictedCluster'] = test_preds.loc[idx, 'PredictedCluster']
        mastertable.loc[idx, 'Above90'] = test_preds.loc[idx, 'Above90']
        mastertable.loc[idx, 'Top70 Perc. Confident'] = test_preds.loc[idx, 'Top70 Perc. Confident']

    if idx in staudt_only_labels.index:
        mastertable.loc[idx, 'Staudt Only Cluster'] = staudt_only_labels.loc[idx, 'cluster']

    if idx in shipp_only_labels.index:
        mastertable.loc[idx, 'Shipp Only Cluster'] = shipp_only_labels.loc[idx, 'cluster']


# Fill in mean/median coverage
mastertable['Mean Target Coverage'] = coverages['coverage_mean']
mastertable['Median Target Coverage'] = coverages['coverage_median']

mastertable = mastertable.sort_values(by='Cluster')
mastertable.to_csv('../../data_tables/tableS1_classifier.tsv', sep='\t', na_rep='NA')