import pandas as pd
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import umap
import random
from sklearn.decomposition import PCA
import seaborn as sns
import src_python.format_data as fd
from skfda import FDataGrid
from skfda.exploratory.stats import geometric_median

seed = 1234
random.seed(seed)
np.random.seed(seed)


def draw_pie(dist,
             xpos,
             ypos,
             size,
             ax=None):
    if ax is None:
        fig, ax = plt.subplots(figsize=(10,8))

    # for incremental pie slices
    cumsum = np.cumsum(dist)
    cumsum = cumsum/ cumsum[-1]
    pie = [0] + cumsum.tolist()

    for r1, r2 in zip(pie[:-1], pie[1:]):
        angles = np.linspace(2 * np.pi * r1, 2 * np.pi * r2)
        x = [0] + np.cos(angles).tolist()
        y = [0] + np.sin(angles).tolist()

        xy = np.column_stack([x, y])

        ax.scatter([xpos], [ypos], marker=xy, s=size)

    return ax


def drawPieMarker(xs, ys, ratios, sizes, colors, ax):
    assert sum(ratios) <= 1, 'sum of ratios needs to be < 1'

    markers = []
    previous = 0
    # calculate the points of the pie pieces
    for color, ratio in zip(colors, ratios):
        this = 2 * np.pi * ratio + previous
        x = [0] + np.cos(np.linspace(previous, this, 10)).tolist() + [0]
        y = [0] + np.sin(np.linspace(previous, this, 10)).tolist() + [0]
        xy = np.column_stack([x, y])
        previous = this
        markers.append({'marker': xy, 's': np.abs(xy).max()**2*np.array(sizes), 'facecolor': color})

    # scatter each of the pie pieces to create pies
    for marker in markers:
        ax.scatter(xs, ys, **marker)


def euclidean_distance(p1, p2):
    return np.linalg.norm(np.array(p1) - np.array(p2))


def detect_outliers(clus_samples, u):
    X = FDataGrid(u.loc[clus_samples])
    median = geometric_median(X)
    median = np.array([median.data_matrix[0][0][0], median.data_matrix[0][1][0]])
    u_clus = u.loc[clus_samples].copy(deep=True)
    u_clus['dis_to_median'] = u_clus.loc[clus_samples].apply(lambda x: euclidean_distance(x, median), axis=1)
    sd = np.std(u_clus['dis_to_median'])
    u_clus['outlier'] = False
    u_clus.loc[u_clus['dis_to_median'] > (3 * sd), 'outlier'] = True
    return u_clus


pallet = {1: sns.color_palette()[4],        # purple, C1
          2: sns.color_palette()[9],        # blue, C2
          3: sns.color_palette()[1],        # orange, C3
          4: sns.color_palette()[2],        # green, C4
          5: sns.color_palette()[3],        # red, C5
          }

preds = '../evaluation_validation_set/confidence_adjusted_tables/NN_reducedV3.3_nfeatures21_pMax0.93344957.tsv'
preds_test = '../evaluation_test_set/NN_reducedV3.3_nfeatures21_testsetEval.tsv'
datafile = '../data_tables/gsm/DLBCL.699.fullGSM.Sep_23_2022.tsv'
training_set = list(pd.read_csv('../data_tables/train_test_sets/TrainingSet_550Subset_May2021.txt', sep='\t', header=None, index_col=0).index)
testing_set = list(pd.read_csv('../data_tables/train_test_sets/TestingSet_149Subset_May2021.txt', sep='\t', header=None, index_col=0).index)
sig_genes_file = '../data_tables/qval_dfs/fisher_exact_5x2.Sep_23_2022.combined.tsv'

df = pd.read_csv(datafile, sep='\t', index_col=0).T
sig_genes = pd.read_csv(sig_genes_file, sep='\t', index_col=0)
sig_genes = sig_genes.loc[sig_genes['q'] <= 0.10].index

# Test+Train UMAP piecharts


def umap_piecharts():
    p_train = pd.read_csv(preds, sep='\t', index_col=0)
    p_test = pd.read_csv(preds_test, sep='\t', index_col=0)

    p_train['cluster'] = p_train.idxmax(axis=1).astype(int) + 1
    p_test['PredictedCluster'] = p_test['PredictedCluster'].astype(int)


    p_train = pd.concat([p_train,
                         pd.DataFrame(p_test[['C1', 'C2', 'C3', 'C4', 'C5', 'PredictedCluster']].values, columns=p_train.columns, index=p_test.index)])
    p_train['Confidence'] = p_train.iloc[:, 0:5].max(axis=1)

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

    # plt.figure(figsize=(10,10), dpi=200)

    fig, ax = plt.subplots(figsize=(5, 5))

    for sample, row in u.iterrows():
        drawPieMarker(xs=row['u1'], ys=row['u2'],
                      ratios=p_train.loc[sample, ['0', '1', '2', '3', '4']], sizes=[100],
                      colors=[pallet[1], pallet[2], pallet[3], pallet[4], pallet[5]],
                      ax=ax)


    plt.savefig('../plots/umap/test_train_umap_piecharts.png', bbox_inches='tight')
    plt.savefig('../plots/umap/test_train_umap_piecharts.pdf', bbox_inches='tight')
    plt.clf()


# Only training set UMAP first
def training_umap_only():
    p_train = pd.read_csv(preds, sep='\t', index_col=0)
    p_test = pd.read_csv(preds_test, sep='\t', index_col=0)
    p_train['cluster'] = p_train.idxmax(axis=1).astype(int) + 1
    p_train['confidence'] = p_train.iloc[:, 0:5].max(axis=1)
    p_train = p_train.loc[p_train['confidence'] > 0.70]
    t_df = df.loc[p_train.index]
    reduced_df = fd.construct_reduced_winning_version(t_df)
    p_train = p_train.loc[reduced_df.index]
    reduced_df.loc[:, 'cluster'] = p_train['cluster']

    to_drop_umap = ['cluster']

    fit = umap.UMAP(n_components=2, n_neighbors=10, min_dist=0.05, metric='euclidean', random_state=seed)
    u = pd.DataFrame(fit.fit_transform(reduced_df.drop(to_drop_umap, axis=1)))
    u.index = reduced_df.index
    u.columns = ['u1', 'u2']

    plt.figure(figsize=(14, 12))

    for clus in range(1, 6):
        clus_samples = reduced_df.loc[reduced_df['cluster'] == clus].index
        plt.scatter(u.loc[clus_samples, 'u1'], u.loc[clus_samples, 'u2'], color=pallet[clus], marker='o', s=40)

    legend_elements = [Line2D([0], [0], marker='_', color=pallet[1], markerfacecolor=pallet[1], label='C1_DLBclass', markersize=12),
                       Line2D([0], [0], marker='_', color=pallet[2], markerfacecolor=pallet[2], label='C2_DLBclass', markersize=12),
                       Line2D([0], [0], marker='_', color=pallet[3], markerfacecolor=pallet[3], label='C3_DLBclass', markersize=12),
                       Line2D([0], [0], marker='_', color=pallet[4], markerfacecolor=pallet[4], label='C4_DLBclass', markersize=12),
                       Line2D([0], [0], marker='_', color=pallet[5], markerfacecolor=pallet[5], label='C5_DLBclass', markersize=12)]

    plt.legend(handles=legend_elements, loc='center left', bbox_to_anchor=(1, 0.5))
    plt.title('Train Set > 0.70 confidence')
    plt.savefig('../plots/umap/train_umap.png', bbox_inches='tight')
    plt.savefig('../plots/umap/train_umap.pdf', bbox_inches='tight')
    plt.clf()

# Test+Train UMAP high confidence spike in

def spikein_umap():
    p_train = pd.read_csv(preds, sep='\t', index_col=0)
    p_test = pd.read_csv(preds_test, sep='\t', index_col=0)

    p_train['cluster'] = p_train.idxmax(axis=1).astype(int) + 1
    p_test['PredictedCluster'] = p_test['PredictedCluster'].astype(int)

    p_all = pd.concat([p_train,
                         pd.DataFrame(p_test[['C1', 'C2', 'C3', 'C4', 'C5', 'PredictedCluster']].values, columns=p_train.columns, index=p_test.index)])
    p_all['Confidence'] = p_all.iloc[:, 0:5].max(axis=1)

    c_cutoff = 0.70
    hc_samples = p_all.loc[p_all['Confidence'] > c_cutoff].index

    reduced_df = fd.construct_reduced_winning_version(df)
    p_all = p_all.loc[reduced_df.index]
    reduced_df.loc[:, 'cluster'] = p_all['cluster']

    reduced_df.loc[:, 'set'] = 'Train'
    reduced_df.loc[testing_set, 'set'] = 'Test'

    # n_spikes = 150
    #spike_in_samples = pd.Series(hc_samples).sample(n=n_spikes, replace=True)
    spike_in_rdf = reduced_df.copy(deep=True)
    spike_in_rdf.index = spike_in_rdf.index + '_spiked'

    reduced_df = pd.concat([reduced_df, spike_in_rdf])

    to_drop_umap = ['set', 'cluster']

    fit = umap.UMAP(n_components=2, n_neighbors=10, min_dist=0.05, metric='euclidean', random_state=seed)
    u = pd.DataFrame(fit.fit_transform(reduced_df.drop(to_drop_umap, axis=1)))
    u.index = reduced_df.index
    u.columns = ['u1', 'u2']

    u = u.loc[p_all.index]
    reduced_df = reduced_df.loc[p_all.index]

    plt.figure(figsize=(14, 12))
    for clus in range(1, 6):
        clus_samples = reduced_df.loc[reduced_df['cluster'] == clus].index

        clus_samples_train = clus_samples[clus_samples.isin(training_set) & clus_samples.isin(hc_samples)]
        clus_samples_train_outlier = clus_samples[clus_samples.isin(training_set) & ~clus_samples.isin(hc_samples)]
        clus_samples_test = clus_samples[clus_samples.isin(testing_set) & clus_samples.isin(hc_samples)]
        clus_samples_test_outlier = clus_samples[clus_samples.isin(testing_set) & ~clus_samples.isin(hc_samples)]

        plt.scatter(u.loc[clus_samples_train, 'u1'], u.loc[clus_samples_train, 'u2'], color=pallet[clus], marker='o', s=50)
        plt.scatter(u.loc[clus_samples_test, 'u1'], u.loc[clus_samples_test, 'u2'], color=pallet[clus], marker='P', s=60)

        plt.scatter(u.loc[clus_samples_train, 'u1'], u.loc[clus_samples_train, 'u2'], s=120, facecolors='none', edgecolors='black')
        plt.scatter(u.loc[clus_samples_test, 'u1'], u.loc[clus_samples_test, 'u2'], s=120, facecolors='none', edgecolors='black', marker='P')

        plt.scatter(u.loc[clus_samples_train_outlier, 'u1'], u.loc[clus_samples_train_outlier, 'u2'], color=pallet[clus], marker='o', s=30)
        plt.scatter(u.loc[clus_samples_test_outlier, 'u1'], u.loc[clus_samples_test_outlier, 'u2'], color=pallet[clus], marker='P', s=40)

    legend_elements = [Line2D([0], [0], marker='_', color=pallet[1], markerfacecolor=pallet[1], label='C1_DLBclass', markersize=12),
                       Line2D([0], [0], marker='_', color=pallet[2], markerfacecolor=pallet[2], label='C2_DLBclass', markersize=12),
                       Line2D([0], [0], marker='_', color=pallet[3], markerfacecolor=pallet[3], label='C3_DLBclass', markersize=12),
                       Line2D([0], [0], marker='_', color=pallet[4], markerfacecolor=pallet[4], label='C4_DLBclass', markersize=12),
                       Line2D([0], [0], marker='_', color=pallet[5], markerfacecolor=pallet[5], label='C5_DLBclass', markersize=12),
                       Line2D([0], [0], marker='o', color='w', markerfacecolor='black', label='Train', markersize=10),
                       Line2D([0], [0], marker='P', color='w', markerfacecolor='black', label='Test', markersize=15)]

    plt.legend(handles=legend_elements, loc='center left', bbox_to_anchor=(1, 0.5))
    plt.savefig('../plots/umap/test_train_umap_spikein_n' + str(699) + '.png', bbox_inches='tight')
    plt.savefig('../plots/umap/test_train_umap_spikein_n' + str(699) + '.pdf', bbox_inches='tight')
    plt.clf()


# Test+Train UMAP rings

def umap_rings():
    p_train = pd.read_csv(preds, sep='\t', index_col=0)
    p_test = pd.read_csv(preds_test, sep='\t', index_col=0)

    p_train['cluster'] = p_train.idxmax(axis=1).astype(int) + 1
    p_test['PredictedCluster'] = p_test['PredictedCluster'].astype(int)

    p_train = pd.concat([p_train,
                         pd.DataFrame(p_test[['C1', 'C2', 'C3', 'C4', 'C5', 'PredictedCluster']].values, columns=p_train.columns, index=p_test.index)])
    p_train['Confidence'] = p_train.iloc[:, 0:5].max(axis=1)

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

    C_CUTOFF = 0.90

    hc_samples = p_train.loc[p_train['Confidence'] > C_CUTOFF].index

    for clus in range(1, 6):
        clus_samples = reduced_df.loc[reduced_df['cluster'] == clus].index

        #u_clus = detect_outliers(clus_samples, u)

        #clus_samples_train = clus_samples[clus_samples.isin(training_set) & ~clus_samples.isin(u_clus.loc[u_clus['outlier']].index)]
        #clus_samples_train_outlier = clus_samples[clus_samples.isin(training_set) & clus_samples.isin(u_clus.loc[u_clus['outlier']].index)]
        #clus_samples_test = clus_samples[clus_samples.isin(testing_set) & ~clus_samples.isin(u_clus.loc[u_clus['outlier']].index)]
        #clus_samples_test_outlier = clus_samples[clus_samples.isin(testing_set) & clus_samples.isin(u_clus.loc[u_clus['outlier']].index)]

        clus_samples_train = clus_samples[clus_samples.isin(training_set) & clus_samples.isin(hc_samples)]
        clus_samples_train_outlier = clus_samples[clus_samples.isin(training_set) & ~clus_samples.isin(hc_samples)]
        clus_samples_test = clus_samples[clus_samples.isin(testing_set) & clus_samples.isin(hc_samples)]
        clus_samples_test_outlier = clus_samples[clus_samples.isin(testing_set) & ~clus_samples.isin(hc_samples)]

        plt.scatter(u.loc[clus_samples_train, 'u1'], u.loc[clus_samples_train, 'u2'], color=pallet[clus], marker='o', s=50)
        plt.scatter(u.loc[clus_samples_test, 'u1'], u.loc[clus_samples_test, 'u2'], color=pallet[clus], marker='P', s=60)

        plt.scatter(u.loc[clus_samples_train, 'u1'], u.loc[clus_samples_train, 'u2'], s=120, facecolors='none', edgecolors='black')
        plt.scatter(u.loc[clus_samples_test, 'u1'], u.loc[clus_samples_test, 'u2'], s=120, facecolors='none', edgecolors='black', marker='P')

        plt.scatter(u.loc[clus_samples_train_outlier, 'u1'], u.loc[clus_samples_train_outlier, 'u2'], color=pallet[clus], marker='o', s=30)
        plt.scatter(u.loc[clus_samples_test_outlier, 'u1'], u.loc[clus_samples_test_outlier, 'u2'], color=pallet[clus], marker='P', s=40)

    legend_elements = [Line2D([0], [0], marker='_', color=pallet[1], markerfacecolor=pallet[1], label='C1_DLBclass', markersize=12),
                       Line2D([0], [0], marker='_', color=pallet[2], markerfacecolor=pallet[2], label='C2_DLBclass', markersize=12),
                       Line2D([0], [0], marker='_', color=pallet[3], markerfacecolor=pallet[3], label='C3_DLBclass', markersize=12),
                       Line2D([0], [0], marker='_', color=pallet[4], markerfacecolor=pallet[4], label='C4_DLBclass', markersize=12),
                       Line2D([0], [0], marker='_', color=pallet[5], markerfacecolor=pallet[5], label='C5_DLBclass', markersize=12),
                       Line2D([0], [0], marker='o', color='w', markerfacecolor='black', label='Train', markersize=10),
                       Line2D([0], [0], marker='P', color='w', markerfacecolor='black', label='Test', markersize=15)]

    plt.legend(handles=legend_elements, loc='center left', bbox_to_anchor=(1, 0.5))
    plt.savefig('../plots/umap/test_train_umap_rings' + str(C_CUTOFF) + 'conf.png', bbox_inches='tight')
    plt.savefig('../plots/umap/test_train_umap_rings' + str(C_CUTOFF) + 'conf.pdf', bbox_inches='tight')
    plt.clf()

# > 70% conf

def umap_hc():
    p_train = pd.read_csv(preds, sep='\t', index_col=0)
    p_test = pd.read_csv(preds_test, sep='\t', index_col=0)

    sig_genes = pd.read_csv(sig_genes_file, sep='\t', index_col=0)
    sig_genes = sig_genes.loc[sig_genes['q'] <= 0.10].index

    p_train['cluster'] = p_train.idxmax(axis=1).astype(int) + 1
    p_test['PredictedCluster'] = p_test['PredictedCluster'].astype(int)

    p_train = pd.concat([p_train,
                         pd.DataFrame(p_test[['C1', 'C2', 'C3', 'C4', 'C5', 'PredictedCluster']].values, columns=p_train.columns, index=p_test.index)])
    #reduced_df = df[sig_genes]

    p_train['confidence'] = p_train.iloc[:, 0:5].max(axis=1)
    p_train = p_train.loc[p_train['confidence'] > 0.70]

    reduced_df = fd.construct_reduced_winning_version(df)
    reduced_df = reduced_df.loc[p_train.index]
    reduced_df.loc[:, 'cluster'] = p_train['cluster']

    reduced_df.loc[:, 'set'] = 'Train'
    reduced_df.loc[reduced_df.index.isin(testing_set), 'set'] = 'Test'

    to_drop_umap = ['set', 'cluster']

    fit = umap.UMAP(n_components=2, n_neighbors=10, min_dist=0.05, metric='euclidean', random_state=seed)
    u = pd.DataFrame(fit.fit_transform(reduced_df.drop(to_drop_umap, axis=1)))
    u.index = reduced_df.index
    u.columns = ['u1', 'u2']

    u.to_csv('../data_tables/umap/umap_70conf.tsv', sep='\t')

    print(p_train['confidence'].min())

    plt.figure(figsize=(14, 12))

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

    plt.title('Train / Test > 0.70 Confidence')
    plt.legend(handles=legend_elements, loc='center left', bbox_to_anchor=(1, 0.5))
    plt.savefig('../plots/umap/test_train_70conf_umap.png', bbox_inches='tight')
    plt.savefig('../plots/umap/test_train_70conf_umap.pdf', bbox_inches='tight')


def umap_custom():
    def distance_dlbcl(x, y, alpha=1):
        c_x = x[-5:]
        y_x = y[-5:]
        f_x = x[0:21]
        f_y = y[0:21]
        dis_f = euclidean_distance(f_x, f_y)
        dis_c = alpha * euclidean_distance(c_x, y_x)
        return dis_f + dis_c

    p_train = pd.read_csv(preds, sep='\t', index_col=0)
    p_test = pd.read_csv(preds_test, sep='\t', index_col=0)

    p_train['cluster'] = p_train.idxmax(axis=1).astype(int) + 1
    p_test['PredictedCluster'] = p_test['PredictedCluster'].astype(int)

    p_all = pd.concat([p_train,
                         pd.DataFrame(p_test[['C1', 'C2', 'C3', 'C4', 'C5', 'PredictedCluster']].values, columns=p_train.columns, index=p_test.index)])
    p_all['Confidence'] = p_all.iloc[:, 0:5].max(axis=1)

    reduced_df = fd.construct_reduced_winning_version(df)
    p_all = p_all.loc[reduced_df.index]
    reduced_df.loc[:, 'cluster'] = p_all['cluster']

    reduced_df.loc[:, 'set'] = 'Train'
    reduced_df.loc[testing_set, 'set'] = 'Test'

    reduced_df[['C1', 'C2', 'C3', 'C4', 'C5']] = p_all[['0', '1', '2', '3', '4']]

    to_drop_umap = ['set', 'cluster']

    fit = umap.UMAP(n_components=2, n_neighbors=10, min_dist=0.05, metric=distance_dlbcl, random_state=seed)
    u = pd.DataFrame(fit.fit_transform(reduced_df.drop(to_drop_umap, axis=1)))
    u.index = reduced_df.index
    u.columns = ['u1', 'u2']

    plt.figure(figsize=(14, 12))

    C_CUTOFF = 0.70

    hc_samples = p_all.loc[p_all['Confidence'] > C_CUTOFF].index

    for clus in range(1, 6):
        clus_samples = reduced_df.loc[reduced_df['cluster'] == clus].index

        clus_samples_train = clus_samples[clus_samples.isin(training_set) & clus_samples.isin(hc_samples)]
        clus_samples_train_outlier = clus_samples[clus_samples.isin(training_set) & ~clus_samples.isin(hc_samples)]
        clus_samples_test = clus_samples[clus_samples.isin(testing_set) & clus_samples.isin(hc_samples)]
        clus_samples_test_outlier = clus_samples[clus_samples.isin(testing_set) & ~clus_samples.isin(hc_samples)]

        plt.scatter(u.loc[clus_samples_train, 'u1'], u.loc[clus_samples_train, 'u2'], color=pallet[clus], marker='o', s=50)
        plt.scatter(u.loc[clus_samples_test, 'u1'], u.loc[clus_samples_test, 'u2'], color=pallet[clus], marker='P', s=60)

        plt.scatter(u.loc[clus_samples_train, 'u1'], u.loc[clus_samples_train, 'u2'], s=120, facecolors='none', edgecolors='black')
        plt.scatter(u.loc[clus_samples_test, 'u1'], u.loc[clus_samples_test, 'u2'], s=120, facecolors='none', edgecolors='black', marker='P')

        plt.scatter(u.loc[clus_samples_train_outlier, 'u1'], u.loc[clus_samples_train_outlier, 'u2'], color=pallet[clus], marker='o', s=30)
        plt.scatter(u.loc[clus_samples_test_outlier, 'u1'], u.loc[clus_samples_test_outlier, 'u2'], color=pallet[clus], marker='P', s=40)

    legend_elements = [Line2D([0], [0], marker='_', color=pallet[1], markerfacecolor=pallet[1], label='C1_DLBclass', markersize=12),
                       Line2D([0], [0], marker='_', color=pallet[2], markerfacecolor=pallet[2], label='C2_DLBclass', markersize=12),
                       Line2D([0], [0], marker='_', color=pallet[3], markerfacecolor=pallet[3], label='C3_DLBclass', markersize=12),
                       Line2D([0], [0], marker='_', color=pallet[4], markerfacecolor=pallet[4], label='C4_DLBclass', markersize=12),
                       Line2D([0], [0], marker='_', color=pallet[5], markerfacecolor=pallet[5], label='C5_DLBclass', markersize=12),
                       Line2D([0], [0], marker='o', color='w', markerfacecolor='black', label='Train', markersize=10),
                       Line2D([0], [0], marker='P', color='w', markerfacecolor='black', label='Test', markersize=15)]

    plt.legend(handles=legend_elements, loc='center left', bbox_to_anchor=(1, 0.5))
    plt.savefig('../plots/umap/test_train_customDis.png', bbox_inches='tight')
    plt.savefig('../plots/umap/test_train_umap_customDis.pdf', bbox_inches='tight')
    plt.clf()


umap_hc()
training_umap_only()