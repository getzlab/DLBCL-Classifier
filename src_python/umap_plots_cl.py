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

targetfile = '../data_tables/confidence_tables/baseline_probabilities.connectivity_based.sensitivity_power2.Sep_23_2022.tsv'
datafile = '../data_tables/gsm/DLBCL.699.fullGSM.Sep_23_2022.tsv'
training_set = list(pd.read_csv('../data_tables/train_test_sets/TrainingSet_550Subset_May2021.txt', sep='\t', header=None, index_col=0).index)
qval_file = '../data_tables/qval_dfs/fisher_exact_5x2.Sep_23_2022.combined.tsv'
ccgd_file = '../data_tables/gsm/GSM.CCGD.updated.Sep_23_2022.tsv'
ccle_file = '../data_tables/gsm/GSM.Depmap.updated.Sep_23_2022.tsv'
cl_mapping_file = '../data_tables/gsm/cell_lines_mappings.tsv'

targets = pd.read_csv(targetfile, sep='\t', index_col=0)
ccgd_preds = pd.read_csv('../evaluation_cell_lines/CCGD_predictions.tsv', sep='\t', index_col=0)
ccle_preds = pd.read_csv('../evaluation_cell_lines/CCLE_predictions.tsv', sep='\t', index_col=0)
ccle_mappings = pd.read_csv(cl_mapping_file, sep='\t', index_col=0)
ccgd_mappings = pd.read_csv(cl_mapping_file, sep='\t', index_col=1)

qval_df = pd.read_csv(qval_file, delimiter='\t', index_col=0)
train_df = pd.read_csv(datafile, sep='\t', index_col=0).T
ccgd_df = pd.read_csv(ccgd_file, sep='\t', index_col=0)
ccle_df = pd.read_csv(ccle_file, sep='\t', index_col=0)

dlbcl_lines = ['DHL2', 'U2932', 'LY4', 'LY7', 'K422', 'BALM3', 'LY19', 'DB', 'WSU_DLBCL2',
               'LY8', 'DHL6', 'WSU_NHL', 'DHL10', 'LY18', 'LY1', 'DHL7', 'DHL16', 'DHL4',
               'PFEIFFER', 'TOLEDO', 'DHL5', 'CTB_1', 'HT', 'DHL8', 'TMD8', 'LY10', 'LY3', 'HBL1']

reduced_train = fd.construct_reduced_winning_version(train_df)
reduced_ccgd = fd.construct_reduced_winning_version(ccgd_df)
reduced_ccle = fd.construct_reduced_winning_version(ccle_df)

targets_train = targets.loc[training_set]
targets_ccgd = ccgd_preds.loc[reduced_ccgd.index]
targets_ccle = ccle_preds.loc[reduced_ccle.index]

reduced_train = reduced_train.loc[training_set]

reduced_train['cluster'] = targets_train['cluster']
reduced_ccgd['cluster'] = targets_ccgd['PredictedCluster']
reduced_ccle['cluster'] = targets_ccle['PredictedCluster']

mappings_ccle = [ccle_mappings.loc[x, 'Ours'] if x in ccle_mappings.index else 'None' for x in targets_ccle.index]
mappings_ccgd = [ccgd_mappings.loc[x, 'DepMap'] if x in ccgd_mappings.index else 'None' for x in targets_ccgd.index]

reduced_ccle['mapped_line'] = mappings_ccle
reduced_train['mapped_line'] = 'None'
reduced_ccgd['mapped_line'] = mappings_ccgd

reduced_train['set'] = 'Train'
reduced_ccgd['set'] = 'CCGD'
reduced_ccle['set'] = 'CCLE'

data_all = pd.concat([reduced_train, reduced_ccgd, reduced_ccle])

pallet = {1: sns.color_palette()[4],        # purple, C1
          2: sns.color_palette()[9],        # blue, C2
          3: sns.color_palette()[1],        # orange, C3
          4: sns.color_palette()[2],        # green, C4
          5: sns.color_palette()[3],        # red, C5
          }

shapes = {'Train': '.',
          'CCLE': 'x',
          'CCGD': 'D'}

to_drop_umap = ['set', 'cluster', 'mapped_line']
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
        if g == 'CCLE' or g == 'CCGD':
            ax.scatter(sub_sub_df['U1'], sub_sub_df['U2'], c=[pallet[clus]], label='C' + str(int(clus)), s=20, marker=shapes[g])
        else:
            ax.scatter(sub_sub_df['U1'], sub_sub_df['U2'], c=[pallet[clus]], label='C' + str(int(clus)), s=15, marker=shapes[g])

legend_elements = [Line2D([0], [0], marker='o', color='w', markerfacecolor=pallet[1], label='C1', markersize=12),
                   Line2D([0], [0], marker='o', color='w', markerfacecolor=pallet[2], label='C2', markersize=12),
                   Line2D([0], [0], marker='o', color='w', markerfacecolor=pallet[3], label='C3', markersize=12),
                   Line2D([0], [0], marker='o', color='w', markerfacecolor=pallet[4], label='C4', markersize=12),
                   Line2D([0], [0], marker='o', color='w', markerfacecolor=pallet[5], label='C5', markersize=12),
                   Line2D([0], [0], marker='X', color='w', markerfacecolor='black', label='CCLE', markersize=12),
                   Line2D([0], [0], marker='D', color='w', markerfacecolor='black', label='CCGD', markersize=12)]

ax.legend(handles=legend_elements)
plt.savefig('../plots/umap/umap_cl_trainset.jpeg')
plt.savefig('../plots/umap/umap_cl_trainset.pdf')
plt.clf()

########################
# DLBCL lines only     #
########################

fig, ax = plt.subplots()
fig.set_size_inches(10, 8)
clusters = set(data_all['cluster'])
grps = set(data_all['set'])
grps.remove('CCLE')
for g in grps:
    sub_df = u.loc[u['set'] == g]
    for clus in clusters:
        sub_sub_df = sub_df.loc[sub_df['cluster'] == clus]
        if g == 'CCGD':
            sub_sub_df = sub_sub_df.loc[sub_sub_df.index.isin(dlbcl_lines)]
            ax.scatter(sub_sub_df['U1'], sub_sub_df['U2'], c=[pallet[clus]], label='C' + str(int(clus)), s=20, marker=shapes[g])
            for idx, row in sub_sub_df.iterrows():
                currtext = idx
                if len(idx) > 8:
                    currtext = idx[0:9]
                ax.annotate(currtext, [sub_sub_df.loc[idx, 'U1'], sub_sub_df.loc[idx, 'U2']])
        else:
            ax.scatter(sub_sub_df['U1'], sub_sub_df['U2'], c=[pallet[clus]], label='C' + str(int(clus)), s=15, marker=shapes[g])

legend_elements = [Line2D([0], [0], marker='o', color='w', markerfacecolor=pallet[1], label='C1', markersize=12),
                   Line2D([0], [0], marker='o', color='w', markerfacecolor=pallet[2], label='C2', markersize=12),
                   Line2D([0], [0], marker='o', color='w', markerfacecolor=pallet[3], label='C3', markersize=12),
                   Line2D([0], [0], marker='o', color='w', markerfacecolor=pallet[4], label='C4', markersize=12),
                   Line2D([0], [0], marker='o', color='w', markerfacecolor=pallet[5], label='C5', markersize=12),
                   Line2D([0], [0], marker='D', color='w', markerfacecolor='black', label='CCGD', markersize=12)]

ax.legend(handles=legend_elements)
plt.savefig('../plots/umap/umap_dlbclLines_trainset.jpeg')
plt.savefig('../plots/umap/umap_dlbclLines_trainset.pdf')
plt.clf()

########################
# High confidence only #
########################
confidences = [0.70, 0.80, 0.90]
ANNOTATE = True
for c_cutoff in confidences:
    ccle_highconf = list(ccle_preds.loc[ccle_preds['Confidence'] > c_cutoff].index)
    ccgd_highconf = list(ccgd_preds.loc[ccgd_preds['Confidence'] > c_cutoff].index)
    trainset_highconf = list(targets_train.loc[targets_train['confidence'] > c_cutoff].index)

    all_high_conf = ccle_highconf + ccgd_highconf + trainset_highconf

    data_highconf = data_all.loc[all_high_conf]

    tmp = data_highconf.loc[data_highconf['set'] == 'CCLE']
    connect_points = []
    for x in tmp.index:
        if x not in ccle_mappings.index:
            continue

        mapped_ccgd = ccle_mappings.loc[x, 'Ours']
        if mapped_ccgd in data_highconf.index:
            connect_points.append((x, mapped_ccgd))

    fit = umap.UMAP(n_components=2, n_neighbors=30, min_dist=0.1, metric='manhattan', random_state=seed)
    u = pd.DataFrame(fit.fit_transform(data_highconf.drop(to_drop_umap, axis=1)))
    u.index = data_highconf.index
    u.columns = ['U1', 'U2']
    u['cluster'] = data_highconf['cluster']
    u['set'] = data_highconf['set']

    fig, ax = plt.subplots()
    fig.set_size_inches(10, 8)
    clusters = set(data_highconf['cluster'])
    grps = set(data_highconf['set'])
    for g in grps:
        sub_df = u.loc[u['set'] == g]
        for clus in clusters:
            sub_sub_df = sub_df.loc[sub_df['cluster'] == clus]
            if g == 'CCLE' or g == 'CCGD':
                ax.scatter(sub_sub_df['U1'], sub_sub_df['U2'], c=[pallet[clus]], label='C' + str(int(clus)), s=35, marker=shapes[g])
                if ANNOTATE:
                    for idx, row in sub_sub_df.iterrows():
                        currtext = idx
                        if len(idx) > 8:
                            currtext = idx[0:9]
                        ax.annotate(currtext, [sub_sub_df.loc[idx, 'U1'], sub_sub_df.loc[idx, 'U2']])
            else:
                ax.scatter(sub_sub_df['U1'], sub_sub_df['U2'], c=[pallet[clus]], label='C' + str(int(clus)), s=15, marker=shapes[g])

    legend_elements = [Line2D([0], [0], marker='o', color='w', markerfacecolor=pallet[1], label='C1', markersize=12),
                       Line2D([0], [0], marker='o', color='w', markerfacecolor=pallet[2], label='C2', markersize=12),
                       Line2D([0], [0], marker='o', color='w', markerfacecolor=pallet[3], label='C3', markersize=12),
                       Line2D([0], [0], marker='o', color='w', markerfacecolor=pallet[4], label='C4', markersize=12),
                       Line2D([0], [0], marker='o', color='w', markerfacecolor=pallet[5], label='C5', markersize=12),
                       Line2D([0], [0], marker='X', color='w', markerfacecolor='black', label='CCLE', markersize=12),
                       Line2D([0], [0], marker='D', color='w', markerfacecolor='black', label='CCGD', markersize=12)]

    for i in range(len(connect_points)):
        curr_pair = connect_points[i]
        ccle_s, ccgd_s = curr_pair[0], curr_pair[1]
        x_ccle, y_ccle = u.loc[ccle_s, 'U1'], u.loc[ccle_s, 'U2']
        x_ccgd, y_ccgd = u.loc[ccgd_s, 'U1'], u.loc[ccgd_s, 'U2']
        ax.plot([x_ccle, x_ccgd], [y_ccle, y_ccgd], c='black')

    ax.legend(handles=legend_elements)
    plt.title('Confidence >= ' + str(c_cutoff))
    plt.savefig('../plots/umap/umap_cl_trainset_conf' + str(c_cutoff) + '.jpeg')
    plt.savefig('../plots/umap/umap_cl_trainset_conf' + str(c_cutoff) + '.pdf')
    plt.clf()

#####################################################
# High confidence + in-common cell lines. Same loop #
#####################################################

    condition = (data_highconf['set'] == 'Train') | (data_highconf['mapped_line'] != 'None')
    data_highconf_commoncl = data_highconf.loc[condition]

    tmp = data_highconf_commoncl.loc[data_highconf_commoncl['set'] == 'CCLE']
    connect_points = [(x, ccle_mappings.loc[x, 'Ours']) for x in tmp.index if ccle_mappings.loc[x, 'Ours'] in data_highconf_commoncl.index]

    fit = umap.UMAP(n_components=2, n_neighbors=30, min_dist=0.1, metric='manhattan', random_state=seed)
    u = pd.DataFrame(fit.fit_transform(data_highconf_commoncl.drop(to_drop_umap, axis=1)))
    u.index = data_highconf_commoncl.index
    u.columns = ['U1', 'U2']
    u['cluster'] = data_highconf_commoncl['cluster']
    u['set'] = data_highconf_commoncl['set']

    fig, ax = plt.subplots()
    fig.set_size_inches(10, 8)
    clusters = set(data_highconf_commoncl['cluster'])
    grps = set(data_highconf_commoncl['set'])
    for g in grps:
        sub_df = u.loc[u['set'] == g]
        for clus in clusters:
            sub_sub_df = sub_df.loc[sub_df['cluster'] == clus]
            if g == 'CCLE' or g == 'CCGD':
                ax.scatter(sub_sub_df['U1'], sub_sub_df['U2'], c=[pallet[clus]], label='C' + str(int(clus)), s=35, marker=shapes[g])
                if ANNOTATE:
                    for idx, row in sub_sub_df.iterrows():
                        currtext = idx
                        if len(idx) > 8:
                            currtext = idx[0:9]
                        ax.annotate(currtext, [sub_sub_df.loc[idx, 'U1'], sub_sub_df.loc[idx, 'U2']])
            else:
                ax.scatter(sub_sub_df['U1'], sub_sub_df['U2'], c=[pallet[clus]], label='C' + str(int(clus)), s=15, marker=shapes[g])

    for i in range(len(connect_points)):
        curr_pair = connect_points[i]
        ccle_s, ccgd_s = curr_pair[0], curr_pair[1]
        x_ccle, y_ccle = u.loc[ccle_s, 'U1'], u.loc[ccle_s, 'U2']
        x_ccgd, y_ccgd = u.loc[ccgd_s, 'U1'], u.loc[ccgd_s, 'U2']
        ax.plot([x_ccle, x_ccgd], [y_ccle, y_ccgd], c='black')

    legend_elements = [Line2D([0], [0], marker='o', color='w', markerfacecolor=pallet[1], label='C1', markersize=12),
                       Line2D([0], [0], marker='o', color='w', markerfacecolor=pallet[2], label='C2', markersize=12),
                       Line2D([0], [0], marker='o', color='w', markerfacecolor=pallet[3], label='C3', markersize=12),
                       Line2D([0], [0], marker='o', color='w', markerfacecolor=pallet[4], label='C4', markersize=12),
                       Line2D([0], [0], marker='o', color='w', markerfacecolor=pallet[5], label='C5', markersize=12),
                       Line2D([0], [0], marker='X', color='w', markerfacecolor='black', label='CCLE', markersize=12),
                       Line2D([0], [0], marker='D', color='w', markerfacecolor='black', label='CCGD', markersize=12)]

    ax.legend(handles=legend_elements)
    plt.title('Confidence >= ' + str(c_cutoff) + ', In-common Cell Lines')
    plt.savefig('../plots/umap/umap_cl_trainset_conf' + str(c_cutoff) + '_commoncl.jpeg')
    plt.savefig('../plots/umap/umap_cl_trainset_conf' + str(c_cutoff) + '_commoncl.pdf')
    plt.clf()



#print(ccgd_mappings.shape)
#data_highconf_mapped =

