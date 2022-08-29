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

targetfile = '../data_tables/confidence_tables/baseline_probabilities.connectivity_based.sensitivity_power2.Aug_17_2022.tsv'
datafile = '../data_tables/gsm/DLBCL_Staudt_Shipp_CL.for_classifier_training.classifier_subset.fix_sv.fix_ploidy.17-Aug-2022.txt'
training_set = list(pd.read_csv('../data_tables/train_test_sets/TrainingSet_550Subset_May2021.txt', sep='\t', header=None, index_col=0).index)
qval_file = '../data_tables/qval_dfs/fisher_exact_5x2_17-Aug-2022.combined.tsv'
lng_file = '../data_tables/gsm/Quality_Longitudinals.01-Sep-2021.txt'

targets = pd.read_csv(targetfile, sep='\t', index_col=0)
lng_preds = pd.read_csv('../evaluation_longitudinals/longitudinals_predictions.tsv', sep='\t', index_col=0)

qval_df = pd.read_csv(qval_file, delimiter='\t', index_col=0)
train_df = pd.read_csv(datafile, sep='\t', index_col=0).T
lng_df = pd.read_csv(lng_file, sep='\t', index_col=0).T

# Fix up longitudinals gsm
lng_df = lng_df.drop('cohort')
lng_df = lng_df.loc[lng_preds.index]
lng_df.columns = ['X' + x if '.AMP' in x or '.DEL' in x else x for x in lng_df.columns]
lng_df.columns = lng_df.columns.str.replace('_', '.').str.replace('-', '.')
lng_df['MYD88'] = lng_df[['MYD88.OTHER', 'MYD88.L265P']].max(axis=1)


reduced_train = fd.construct_reduced_winning_version(train_df)
reduced_lng = fd.construct_reduced_winning_version(lng_df)

targets_train = targets.loc[training_set]
preds_lng = lng_preds.loc[reduced_lng.index]

reduced_train = reduced_train.loc[training_set]

reduced_train['cluster'] = targets_train['cluster']
reduced_lng['cluster'] = preds_lng['PredictedCluster']

next_sample_map = {'DFCIDL001_DT': 'DFCIDL001_R1',
                   'DFCIDL003_DT': 'DFCIDL003_R1',
                   'DFCIDL004_DT': 'DFCIDL004_R1',
                   'DFCIDL007_DT': 'DFCIDL007_R1',
                   'DFCIDL008_DT': 'DFCIDL008_R1', 'DFCIDL008_R1': 'DFCIDL008_R2',
                   'DFCIDL009_DT': 'DFCIDL009_R1', 'DFCIDL009_R1': 'DFCIDL009_R2',
                   'DFCIDL010_DT': 'DFCIDL010_R1',
                   'DLBCL_c_D_pair10_R1': 'DLBCL_c_D_pair10_R2',
                   'DLBCL_c_D_pair13': 'DLBCL_c_D_pair13_R1',
                   'DLBCL_c_D_pair20': 'DLBCL_c_D_pair20_R1', 'DLBCL_c_D_pair20_R1': 'DLBCL_c_D_pair20_R2',
                   'DLBCL_c_D_pair22': 'DLBCL_c_D_pair22_R1',
                   'DLBCL_c_D_pair23': 'DLBCL_c_D_pair23_R1',
                   'DLBCL_c_D_pair2': 'DLBCL_c_D_pair2_R1',
                   'DLBCL_c_D_pair5': 'DLBCL_c_D_pair5_R1',
                   'DLBCL_c_D_pair8': 'DLBCL_c_D_pair8_R1'}


reduced_train['set'] = 'Train Set'
reduced_lng['set'] = 'Relapse'

data_all = pd.concat([reduced_train, reduced_lng])
data_all = data_all[~data_all.index.duplicated(keep='last')]

pallet = {1: sns.color_palette()[4],        # purple, C1
          2: sns.color_palette()[9],        # blue, C2
          3: sns.color_palette()[1],        # orange, C3
          4: sns.color_palette()[2],        # green, C4
          5: sns.color_palette()[3],        # red, C5
          }

shapes = {'Train Set': '.',
          'Relapse': 'x',
          'Relapse High Confidence': 's'}

to_drop_umap = ['set', 'cluster']
################################
# All samples reduced features #
################################

fit = umap.UMAP(n_components=2, n_neighbors=30, min_dist=0.1, metric='manhattan', random_state=seed)
u = pd.DataFrame(fit.fit_transform(data_all.drop(to_drop_umap, axis=1)))
u.index = data_all.index
u.columns = ['U1', 'U2']
u['cluster'] = data_all['cluster']
u['set'] = data_all['set']

fig, ax = plt.subplots()
fig.set_size_inches(10, 8)
clusters = set(data_all['cluster'])
grps = set(data_all['set'])
for g in grps:
    sub_df = u.loc[u['set'] == g]
    for clus in clusters:
        sub_sub_df = sub_df.loc[sub_df['cluster'] == clus]
        if 'Relapse' in g:
            ax.scatter(sub_sub_df['U1'], sub_sub_df['U2'], c=[pallet[clus]], label='C' + str(int(clus)), s=60, marker=shapes[g])
            for idx, row in sub_sub_df.iterrows():
                currtext = idx
                ax.annotate(currtext, [sub_sub_df.loc[idx, 'U1'], sub_sub_df.loc[idx, 'U2']], fontsize=6)
        else:
            ax.scatter(sub_sub_df['U1'], sub_sub_df['U2'], c=[pallet[clus]], label='C' + str(int(clus)), s=15, marker=shapes[g])

for s in next_sample_map:
    next_s = next_sample_map[s]
    s_x, s_y = u.loc[s, 'U1'], u.loc[s, 'U2']
    ns_x, ns_y = u.loc[next_s, 'U1'], u.loc[next_s, 'U2']
    ax.arrow(x=s_x, y=s_y,
             dx=(ns_x - s_x) / 2, dy=(ns_y - s_y) / 2,
             width=0.00, length_includes_head=False, head_width=0.06, head_length=0.06)
    ax.plot([s_x, ns_x], [s_y, ns_y], c='black')


legend_elements = [Line2D([0], [0], marker='o', color='w', markerfacecolor=pallet[1], label='C1', markersize=12),
                   Line2D([0], [0], marker='o', color='w', markerfacecolor=pallet[2], label='C2', markersize=12),
                   Line2D([0], [0], marker='o', color='w', markerfacecolor=pallet[3], label='C3', markersize=12),
                   Line2D([0], [0], marker='o', color='w', markerfacecolor=pallet[4], label='C4', markersize=12),
                   Line2D([0], [0], marker='o', color='w', markerfacecolor=pallet[5], label='C5', markersize=12),
                   Line2D([0], [0], marker='X', color='w', markerfacecolor='black', label='Relapse Sets', markersize=12)]

ax.legend(handles=legend_elements, bbox_to_anchor=(1, 1))
plt.savefig('../plots/umap/umap_longitudinals_trainset.jpeg')
plt.savefig('../plots/umap/umap_longitudinals_trainset.pdf')
plt.clf()

# ##########################
# # High confidence symbol #
# ##########################

hc_train = list(targets_train.loc[targets_train['confidence'] >= 0.90].index)
hc_lngs = list(preds_lng.loc[preds_lng['Confidence'] >= 0.90].index)
hc_all = hc_train + hc_lngs

samples = sorted(list(set(training_set + list(next_sample_map.keys()) + list(next_sample_map.values()))))
data_paired = data_all.loc[samples]
hc_lngs = [x for x in hc_lngs if x in data_paired.index]
data_paired.loc[hc_lngs, 'set'] = 'Relapse High Confidence'

data_paired = data_paired.loc[data_paired.drop(to_drop_umap, axis=1).astype(bool).sum(axis=1) > 1]

fit = umap.UMAP(n_components=2, n_neighbors=30, min_dist=0.1, metric='manhattan', random_state=seed)
u = pd.DataFrame(fit.fit_transform(data_paired.drop(to_drop_umap, axis=1)))
u.index = data_paired.index
u.columns = ['U1', 'U2']
u['cluster'] = data_paired['cluster']
u['set'] = data_paired['set']

fig, ax = plt.subplots()
fig.set_size_inches(10, 8)
clusters = set(data_paired['cluster'])
grps = set(data_paired['set'])
for g in grps:
    sub_df = u.loc[u['set'] == g]
    for clus in clusters:
        sub_sub_df = sub_df.loc[sub_df['cluster'] == clus]
        if 'Relapse' in g:
            ax.scatter(sub_sub_df['U1'], sub_sub_df['U2'], c=[pallet[clus]], label='C' + str(int(clus)), s=75, marker=shapes[g])
            for idx, row in sub_sub_df.iterrows():
                currtext = idx
                if 'R1' in currtext and currtext != 'DLBCL_c_D_pair10_R1':
                    currtext = 'R1'
                if 'R2' in currtext:
                    currtext = 'R2'
                ax.annotate(currtext, [sub_sub_df.loc[idx, 'U1'] + 0.02, sub_sub_df.loc[idx, 'U2'] + 0.02], fontsize=6)
        else:
            ax.scatter(sub_sub_df['U1'], sub_sub_df['U2'], c=[pallet[clus]], label='C' + str(int(clus)), s=15, marker=shapes[g])

for s in next_sample_map:
    if s not in u.index:
        continue

    next_s = next_sample_map[s]
    s_x, s_y = u.loc[s, 'U1'], u.loc[s, 'U2']
    ns_x, ns_y = u.loc[next_s, 'U1'], u.loc[next_s, 'U2']
    ax.arrow(x=s_x, y=s_y,
             dx=(ns_x - s_x) / 2, dy=(ns_y - s_y) / 2,
             width=0.00, length_includes_head=False, head_width=0.06, head_length=0.06)
    ax.plot([s_x, ns_x], [s_y, ns_y], c='black')


legend_elements = [Line2D([0], [0], marker='o', color='w', markerfacecolor=pallet[1], label='C1', markersize=12),
                   Line2D([0], [0], marker='o', color='w', markerfacecolor=pallet[2], label='C2', markersize=12),
                   Line2D([0], [0], marker='o', color='w', markerfacecolor=pallet[3], label='C3', markersize=12),
                   Line2D([0], [0], marker='o', color='w', markerfacecolor=pallet[4], label='C4', markersize=12),
                   Line2D([0], [0], marker='o', color='w', markerfacecolor=pallet[5], label='C5', markersize=12),
                   Line2D([0], [0], marker='X', color='w', markerfacecolor='black', label='Relapse Low Confidence', markersize=12),
                   Line2D([0], [0], marker='s', color='w', markerfacecolor='black', label='Relapse High Confidence', markersize=12)]

ax.legend(handles=legend_elements, bbox_to_anchor=(1, 1))
plt.savefig('../plots/umap/umap_longitudinals_hcmarker_trainset.jpeg', bbox_inches='tight')
plt.savefig('../plots/umap/umap_longitudinals_hcmarker_trainset.pdf', bbox_inches='tight')
plt.clf()
