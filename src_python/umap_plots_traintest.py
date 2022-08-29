import pandas as pd
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import umap
import random
from sklearn.decomposition import PCA
import seaborn as sns
import format_data as fd

seed = 1234
random.seed(seed)
np.random.seed(seed)

pallet = {1: sns.color_palette()[4],        # purple, C1
          2: sns.color_palette()[9],        # blue, C2
          3: sns.color_palette()[1],        # orange, C3
          4: sns.color_palette()[2],        # green, C4
          5: sns.color_palette()[3],        # red, C5
          }

preds = '../evaluation_validation_set/confidence_adjusted_tables/NN_reducedV3.1_qval0.05_removeN5_nfeatures21_pMax0.9393394.tsv'
preds_test = '../evaluation_test_set/NN_reducedV3.1_qval0.05_removeN5_nfeatures21_testsetEval.tsv'
datafile = '../data_tables/gsm/DLBCL_Staudt_Shipp_CL.for_classifier_training.classifier_subset.fix_sv.fix_ploidy.28-Mar-2022.txt'
training_set = list(pd.read_csv('../data_tables/train_test_sets/TrainingSet_550Subset_May2021.txt', sep='\t', header=None, index_col=0).index)
testing_set = list(pd.read_csv('../data_tables/train_test_sets/TestingSet_149Subset_May2021.txt', sep='\t', header=None, index_col=0).index)
sig_genes_file = '../data_tables/qval_dfs/fisher_exact_5x2_28-Mar-2022.combined.tsv'

df = pd.read_csv(datafile, sep='\t', index_col=0).T
p_train = pd.read_csv(preds, sep='\t', index_col=0)
p_test = pd.read_csv(preds_test, sep='\t', index_col=0)

sig_genes = pd.read_csv(sig_genes_file, sep='\t', index_col=0)
sig_genes = sig_genes.loc[sig_genes['q'] <= 0.05].index

p_train['cluster'] = p_train.idxmax(axis=1).astype(int) + 1
p_test['PredictedCluster'] = p_test['PredictedCluster'].astype(int)

p_train = pd.concat([p_train,
                     pd.DataFrame(p_test[['C1', 'C2', 'C3', 'C4', 'C5', 'PredictedCluster']].values, columns=p_train.columns, index=p_test.index)])
#reduced_df = df[sig_genes]

reduced_df = fd.construct_reduced_winning_version(df)
p_train = p_train.loc[reduced_df.index]
reduced_df.loc[:, 'cluster'] = p_train['cluster']

reduced_df.loc[:, 'set'] = 'Train'
reduced_df.loc[testing_set, 'set'] = 'Test'

to_drop_umap = ['set', 'cluster']

fit = umap.UMAP(n_components=2, n_neighbors=10, min_dist=0.05, metric='euclidean', random_state=seed)
u = pd.DataFrame(fit.fit_transform(reduced_df.drop(to_drop_umap, axis=1)))
u.index = reduced_df.index
u.columns = ['u1', 'u2']

plt.figure(figsize=(14, 12))

print(reduced_df)
for clus in range(1, 6):
    clus_samples = reduced_df.loc[reduced_df['cluster'] == clus].index
    clus_samples_train = clus_samples[clus_samples.isin(training_set)]
    clus_samples_test = clus_samples[clus_samples.isin(testing_set)]

    plt.scatter(u.loc[clus_samples_train, 'u1'], u.loc[clus_samples_train, 'u2'], color=pallet[clus], marker='o', s=40)
    plt.scatter(u.loc[clus_samples_test, 'u1'], u.loc[clus_samples_test, 'u2'], color=pallet[clus], marker='x', s=40)

legend_elements = [Line2D([0], [0], marker='_', color=pallet[1], markerfacecolor=pallet[1], label='C1_DLBclass', markersize=12),
                   Line2D([0], [0], marker='_', color=pallet[2], markerfacecolor=pallet[2], label='C2_DLBclass', markersize=12),
                   Line2D([0], [0], marker='_', color=pallet[3], markerfacecolor=pallet[3], label='C3_DLBclass', markersize=12),
                   Line2D([0], [0], marker='_', color=pallet[4], markerfacecolor=pallet[4], label='C4_DLBclass', markersize=12),
                   Line2D([0], [0], marker='_', color=pallet[5], markerfacecolor=pallet[5], label='C5_DLBclass', markersize=12),
                   Line2D([0], [0], marker='o', color='w', markerfacecolor='black', label='Train', markersize=10),
                   Line2D([0], [0], marker='X', color='w', markerfacecolor='black', label='Test', markersize=15)]

plt.legend(handles=legend_elements, loc='center left', bbox_to_anchor=(1, 0.5))
plt.savefig('../plots/umap/test_train_umap.png', bbox_inches='tight')
plt.savefig('../plots/umap/test_train_umap.pdf', bbox_inches='tight')
