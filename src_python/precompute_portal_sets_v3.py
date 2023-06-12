import pandas as pd
import numpy as np
import umap
import random
import seaborn as sns
import format_data as fd
import plotly.express as px
import json
import plotly.graph_objects as go
from matplotlib import pyplot as plt
from PIL import ImageColor

seed = 1234
random.seed(seed)
np.random.seed(seed)

palette = sns.color_palette().as_hex()
palette = {1: palette[4],  # purple, C1
           2: palette[9],  # blue, C2
           3: palette[1],  # orange, C3
           4: palette[2],  # green, C4
           5: palette[3],  # red, C5
           }


def make_icomut_json(heatmap_gsm, qvals, sample_clusters):
    heatmap_gsm = heatmap_gsm.astype(float).astype(int)

    amps = [x for x in heatmap_gsm.index if '.AMP' in x]
    dels = [x for x in heatmap_gsm.index if '.DEL' in x]

    # add qvals at the end
    qvals = qvals.loc[heatmap_gsm.index]

    for idx in heatmap_gsm.index:
        if idx in amps:
            heatmap_gsm.loc[idx] = heatmap_gsm.loc[idx].map({0: 0, 1: 4, 2: 5})

        if idx in dels:
            heatmap_gsm.loc[idx] = heatmap_gsm.loc[idx].map({0: 0, 1: 6, 2: 7})

    heatmap_gsm.index = 'hmap_' + heatmap_gsm.index
    heatmap_gsm['all_q'] = qvals.values
    heatmap_gsm.loc['cluster'] = sample_clusters.loc['cluster']
    heatmap_gsm.loc['confidence'] = sample_clusters.loc['confidence'].astype(float).round(3)
    heatmap_gsm.loc['set'] = sample_clusters.loc['set']

    heatmap_gsm = heatmap_gsm.fillna('NA')
    heatmap_dict = heatmap_gsm.to_dict()
    icm_json = {'data': heatmap_dict, 'samples': list(heatmap_gsm.columns), 'sets': list(heatmap_gsm.index)}
    return icm_json


def get_umap_distances(umap_df, centroids):
    umap_df['distance_from_centroid'] = np.nan
    umap_df['num_std'] = np.nan

    for g in umap_df.groupby('cluster'):
        distances = np.sqrt((umap_df.loc[umap_df['cluster'] == g[0], 'U1'] - centroids[g[0]][0]) ** 2 +
                            (umap_df.loc[umap_df['cluster'] == g[0], 'U2'] - centroids[g[0]][1]) ** 2)
        std_dis = np.std(distances)
        num_stds = distances / std_dis
        umap_df.loc[distances.index, 'distance_from_centroid'] = distances
        umap_df.loc[distances.index, 'num_std'] = num_stds

    return umap_df


datafile = "../data_tables/gsm/DLBCL.699.fullGSM.Sep_23_2022.tsv"
qval_file = '../data_tables/qval_dfs/fisher_exact_5x2.Sep_23_2022.combined.tsv'
targetfile = "../data_tables/confidence_tables/baseline_probabilities.connectivity_based.sensitivity_power2.Sep_23_2022.tsv"
setfile = '../data_tables/sample_sets/sample_inclusion_table.tsv'
train_preds = '../evaluation_validation_set/confidence_adjusted_tables/NN_reducedV3.4_removeN5_nfeatures21_pMax0.93856484.tsv'
test_preds = '../evaluation_test_set/NN_reducedV3.4_removeN5_nfeatures21_testsetEval.tsv'

test_samples = list(pd.read_csv("../data_tables/train_test_sets/TestingSet_149Subset_May2021.txt",
                                    sep='\t', header=None, index_col=0).index)
train_samples = list(pd.read_csv("../data_tables/train_test_sets/TrainingSet_550Subset_May2021.txt",
                                    sep='\t', header=None, index_col=0).index)

train_confs = pd.read_csv(train_preds, sep='\t', index_col=0)
test_confs = pd.read_csv(test_preds, sep='\t', index_col=0)

train_confs.columns = ['C1', 'C2', 'C3', 'C4', 'C5']
train_confs['Confidence'] = train_confs.max(axis=1)
train_confs['PredictedCluster'] = train_confs.iloc[:, 0:5].idxmax(axis=1)

test_confs = test_confs.drop(['Correctness', 'TrueCluster'], axis=1)
test_confs['PredictedCluster'] = test_confs['PredictedCluster'].map({1.0: 'C1', 2.0: 'C2', 3.0: 'C3', 4.0: 'C4', 5.0: 'C5'})


targets = pd.read_csv(targetfile, sep='\t', index_col=0)
targets['cluster'] = targets['cluster'].map({1: 'C1', 2: 'C2', 3: 'C3', 4: 'C4', 5: 'C5'})

c1_s = targets.loc[targets['cluster'] == 'C1'].sort_values('confidence')
c2_s = targets.loc[targets['cluster'] == 'C2'].sort_values('confidence')
c3_s = targets.loc[targets['cluster'] == 'C3'].sort_values('confidence')
c4_s = targets.loc[targets['cluster'] == 'C4'].sort_values('confidence')
c5_s = targets.loc[targets['cluster'] == 'C5'].sort_values('confidence')

sample_sets = pd.read_csv(setfile, sep='\t', index_col=0)
sample_sets = sample_sets.loc[sample_sets['included_in_clustering']]

shipp_s = sample_sets.loc[sample_sets['cohort'] == 'Shipp'].index
staudt_s = sample_sets.loc[sample_sets['cohort'] == 'Staudt'].index
tcga_s = sample_sets.loc[sample_sets['cohort'] == 'TCGA'].index

labels_test = targets.loc[test_samples]
labels_test = labels_test.sort_values(by=['cluster', 'confidence'], ascending=[True, True])

labels_train = targets.loc[train_samples]
labels_train = labels_train.sort_values(by=['cluster', 'confidence'], ascending=[True, True])

df = pd.read_csv(datafile, sep='\t', index_col=0).T

qval_df = pd.read_csv(qval_file, sep='\t', index_col=0)

reduced_features = qval_df.loc[qval_df['q'] <= 0.10].index
df = df[reduced_features]

df_shipp = df.loc[shipp_s].T
df_staudt = df.loc[staudt_s].T
df_test = df.loc[labels_test.index].T
df_train = df.loc[labels_train.index].T

qv = qval_df.loc[df.columns, 'q']

df_shipp.loc['set'] = 'Chapuy'
df_staudt.loc['set'] = 'Schmitz'
df_test.loc['set'] = 'Test'
df_train.loc['set'] = 'Train'

df_shipp.loc['cluster'] = targets.loc[df_shipp.columns, 'cluster']
df_shipp.loc['confidence'] = targets.loc[df_shipp.columns, 'confidence']

df_staudt.loc['cluster'] = targets.loc[df_staudt.columns, 'cluster']
df_staudt.loc['confidence'] = targets.loc[df_staudt.columns, 'confidence']

df_test.loc['cluster'] = targets.loc[df_test.columns, 'cluster']
df_test.loc['confidence'] = targets.loc[df_test.columns, 'confidence']

df_train.loc['cluster'] = targets.loc[df_train.columns, 'cluster']
df_train.loc['confidence'] = targets.loc[df_train.columns, 'confidence']

df_full = pd.concat([df_test, df_train], axis=1)
df_full.loc['set'] = 'Full Set'

shipp_preds = pd.concat([train_confs.loc[train_confs.index.isin(shipp_s)], test_confs.loc[test_confs.index.isin(shipp_s)]])
staudt_preds = pd.concat([train_confs.loc[train_confs.index.isin(staudt_s)], test_confs.loc[test_confs.index.isin(staudt_s)]])
full_preds = pd.concat([train_confs, test_confs])

# .sort_values(['cluster', 'set', 'confidence'], ascending=[True, True, True], axis=1)

##############
# Added sets #
##############

ccle_full_gsm = pd.read_csv('../data_tables/gsm/GSM.Depmap.updated.Sep_23_2022.tsv', sep='\t', index_col=0)
ccle_preds = pd.read_csv('../evaluation_cell_lines/CCLE_predictions.tsv', sep='\t', index_col=0)
ccle_reduced_gsm = pd.read_csv('../evaluation_cell_lines/CCLE_reducedV3.4_removeN5_gsm.tsv', sep='\t', index_col=0)
ccgd_full_gsm = pd.read_csv('../data_tables/gsm/GSM.CCGD.updated.Sep_23_2022.tsv', sep='\t', index_col=0)
ccgd_preds = pd.read_csv('../evaluation_cell_lines/CCGD_predictions.tsv', sep='\t', index_col=0)
ccgd_reduced_gsm = pd.read_csv('../evaluation_cell_lines/CCGD_reducedV3.4_removeN5_gsm.tsv', sep='\t', index_col=0)
false_neg_50_gsm = '../random_dropout_experiment/experiment_gsms/tmpRandomDroppedGSM_fullfeatures_step0.5.txt'
false_pos_50_gsm = '../random_add_in_experiment/experiment_gsms/RandomAddInGSM_fullfeatures_step0.5.txt'
false_neg_50_preds = '../random_dropout_experiment/preds/dropout_preds_0.5.txt'
false_pos_50_preds = '../random_add_in_experiment/preds/addin_preds_0.5.txt'

dropout_preds = pd.read_csv(false_neg_50_preds, sep='\t', index_col=0)
addin_preds = pd.read_csv(false_pos_50_preds, sep='\t', index_col=0)
dropout_preds['PredictedCluster'] = dropout_preds['Predicted Cluster'] + 1
addin_preds['PredictedCluster'] = addin_preds['Predicted Cluster'] + 1

dropout_preds['PredictedCluster'] = dropout_preds['PredictedCluster'].map({1: 'C1', 2: 'C2', 3: 'C3', 4: 'C4', 5: 'C5'})
addin_preds['PredictedCluster'] = addin_preds['PredictedCluster'].map({1: 'C1', 2: 'C2', 3: 'C3', 4: 'C4', 5: 'C5'})

#dropout_preds['cluster'] = dropout_preds['PredictedCluster']
#addin_preds['cluster'] = addin_preds['PredictedCluster']

df_dropped = pd.read_csv(false_neg_50_gsm, sep='\t', index_col=0).T
df_added = pd.read_csv(false_pos_50_gsm, sep='\t', index_col=0).T

ccgd_preds['PredictedCluster'] = ccgd_preds['PredictedCluster'].map({1.0: 'C1', 2.0: 'C2', 3.0: 'C3', 4.0: 'C4', 5.0: 'C5'})
ccle_preds['PredictedCluster'] = ccle_preds['PredictedCluster'].map({1.0: 'C1', 2.0: 'C2', 3.0: 'C3', 4.0: 'C4', 5.0: 'C5'})

df_ccle = ccle_full_gsm.loc[reduced_features, ccle_preds.index]
df_ccgd = ccgd_full_gsm.loc[reduced_features, ccgd_preds.index]

df_ccle.loc['set'] = 'cell Lines CCLE'
df_ccle.loc['cluster'] = ccle_preds.loc[df_ccle.columns, 'PredictedCluster']
df_ccle.loc['confidence'] = ccle_preds.loc[df_ccle.columns, 'Confidence']

df_ccgd.loc['set'] = 'cell Lines CCGD'
df_ccgd.loc['cluster'] = ccgd_preds.loc[df_ccgd.columns, 'PredictedCluster']
df_ccgd.loc['confidence'] = ccgd_preds.loc[df_ccgd.columns, 'Confidence']

df_added = df_added[reduced_features]
df_dropped = df_dropped[reduced_features]

sample_size = 20
df_added_sub = df_added.sample(sample_size).T
df_dropped_sub = df_dropped.sample(sample_size).T

df_added_sub.loc['set'] = 'fp sim 50%'
df_added_sub.loc['cluster'] = addin_preds.loc[df_added_sub.columns, 'PredictedCluster']
df_added_sub.loc['confidence'] = addin_preds.loc[df_added_sub.columns, 'Confidence']
addin_preds = addin_preds.loc[df_added_sub.columns]
df_added_sub.columns = df_added_sub.columns + '_addin'
addin_preds.index = addin_preds.index + '_addin'

df_dropped_sub.loc['set'] = 'fn sim 50%'
df_dropped_sub.loc['cluster'] = dropout_preds.loc[df_dropped_sub.columns, 'PredictedCluster']
df_dropped_sub.loc['confidence'] = dropout_preds.loc[df_dropped_sub.columns, 'Confidence']
dropout_preds = dropout_preds.loc[df_dropped_sub.columns]
df_dropped_sub.columns = df_dropped_sub.columns + '_dropout'
dropout_preds.index = dropout_preds.index + '_dropout'

# m x n precomputes

base_sets = [(df_shipp, 'Chapuy', shipp_preds), (df_staudt, 'Schmitz', staudt_preds), (df_full, 'Full', full_preds),
             (df_test, 'Test', test_confs), (df_train, 'Train', train_confs)]

df_test_2 = df_test.copy(deep=True)
df_test_2.columns = df_test_2.columns + '_2'
test_confs_2 = test_confs.copy(deep=True)
test_confs_2.index = test_confs_2.index + '_2'

df_staudt_2 = df_staudt.copy(deep=True)
df_staudt_2.columns = df_staudt_2.columns + '_2'
staudt_preds_2 = staudt_preds.copy(deep=True)
staudt_preds_2.index = staudt_preds_2.index + '_2'

added_sets = [(df_ccle, 'CCLE', ccle_preds), (df_ccgd, 'CCGD', ccgd_preds),
              (df_added_sub, 'Add-In', addin_preds), (df_dropped_sub, 'Drop-out', dropout_preds),
              (df_test_2, 'Test2', test_confs_2), (df_staudt_2, 'Schmitz2', staudt_preds_2)]

for bset, name, bpreds in base_sets:
    for aset, name2, apreds in added_sets:

        print('------')
        print(name, name2)
        print('------')

        concat_df = pd.concat([bset, aset], axis=1)
        concat_df = concat_df.sort_values(by=['cluster', 'set', 'confidence'], ascending=[True, True, True], axis=1)
        concat_icm_json = make_icomut_json(concat_df.drop(['cluster', 'set', 'confidence']), qv, concat_df.loc[['cluster', 'set', 'confidence']])

        apreds['set'] = name2
        bpreds['set'] = name

        reduced_concat_df = fd.construct_reduced_winning_version(concat_df.drop(['set', 'cluster', 'confidence'], axis=0), add_missing_features=True)
        pred_df = apreds[['C1', 'C2', 'C3', 'C4', 'C5', 'Confidence', 'PredictedCluster', 'set']]
        pred_df = pd.concat([pred_df, bpreds])

        pred_df = pred_df.sort_values(by=['set', 'PredictedCluster', 'Confidence'], ascending=[True, True, True])

        n_comp = 2
        fit = umap.UMAP(n_components=n_comp, n_neighbors=15, min_dist=0.10, metric='euclidean', random_state=seed)
        u = pd.DataFrame(fit.fit_transform(reduced_concat_df))
        u.index = reduced_concat_df.index
        u.columns = ['U' + str(i) for i in range(1, n_comp + 1)]

        u['sample'] = u.index
        u['cluster'] = concat_df.loc['cluster']

        centroids = {g[0]: (g[1]['U1'].mean(), g[1]['U2'].mean()) for g in u.groupby('cluster')}
        u = get_umap_distances(u, centroids)

        pred_df.loc[:, 'U1'] = u.loc[pred_df.index, 'U1']
        pred_df.loc[:, 'U2'] = u.loc[pred_df.index, 'U2']
        pred_df.index.name = 'sample'
        pred_df.to_csv('/Users/twood/Desktop/dlbcl-classifier-portal/files/classifications/' + name + '_' + name2 + '_predictions.tsv', sep='\t', float_format='%.4f')

        fn = '/Users/twood/Desktop/dlbcl-classifier-portal/files/data/dlbcl_' + name + '_' + name2 + '_data.json'
        fn2 = '../portal_icomut_jsons/dlbcl_' + name + '_' + name2 + '_data.json'
        with open(fn, 'w+') as f:
            json.dump(concat_icm_json, f, indent=4)
        with open(fn2, 'w+') as f:
            json.dump(concat_icm_json, f, indent=4)

        with open('/Users/twood/Desktop/dlbcl-classifier-portal/files/config/config.json', 'r') as f:
            cfg_file = json.load(f)

        cfg_file['highlight'] = {'color': '#fdffde', 'samples': list(aset.columns)}
        del cfg_file['samplesToHighlight']
        del cfg_file['highlightColor']
        cfg_file['panels'][0]['label'] = "Cluster"

        set_panel = {
                        "type": "strip",
                        "id": "Set",
                        "label": "Set",
                        "dataPrefix": "set",
                        "row": 1,
                        "colors": [
                            "#000000",
                            "#f5f242"
                        ],
                        "plots": {
                            "center": {
                                "yAxisTitle": "Set"
                            }
                        }
                    }
        confidence_panel = {
                                "type": "continuousstrip",
                                "id": "Confidence",
                                "label": "Confidence",
                                "dataPrefix": "confidence",
                                "row": 2,
                                "colors": [
                                    "#ff1900",
                                    "#00fa08"
                                ],
                                "plots": {
                                    "center": {
                                        "yAxisTitle": "Confidence",
                                        "dataOptions": {"colorRange": {"min": 0.20, "max": 1.0}}
                                    }
                                }
                            }

        cfg_file['panels'].insert(1, set_panel)
        cfg_file['panels'].insert(2, confidence_panel)
        cfg_file['panels'][3]['row'] = 3

        for i in range(0, len(cfg_file['panels'])):
            cfg_file['panels'][i]['displayLegend'] = True

        cfg_file['panels'][0]['displayLegend'] = False

        cfg_file['autoAddPadding'] = True

        cfg_file['initialSort'] = [{"panel": "Cluster", "order": 1},
                                   {"panel": "Set", "order": 1},
                                   {"panel": "Confidence", "order": 1}]

        fn = '/Users/twood/Desktop/dlbcl-classifier-portal/files/config/config_' + name + '_' + name2 + '.json'
        fn2 = '../portal_config_jsons/config_' + name + '_' + name2 + '.json'
        with open(fn, 'w+') as f:
            json.dump(cfg_file, f, indent=4)
        with open(fn2, 'w+') as f:
            json.dump(cfg_file, f, indent=4)


        #########
        # UMAPs #
        #########

        samples = aset.columns

        fig = go.Figure()

        for clus in ['C1', 'C2', 'C3', 'C4', 'C5']:
            selection_cond_bset = ~u.index.isin(samples) & (u['cluster'] == clus)
            selection_cond_aset = u.index.isin(samples) & (u['cluster'] == clus)

            fig.add_trace(
                go.Scatter(mode='markers',
                           x=u.loc[selection_cond_bset, 'U1'],
                           y=u.loc[selection_cond_bset, 'U2'],
                           marker=dict(
                               color=u.loc[selection_cond_bset, 'cluster'].map(
                                   {'C1': palette[1], 'C2': palette[2], 'C3': palette[3], 'C4': palette[4], 'C5': palette[5]}
                               ),
                               size=10,
                               opacity=0.3,
                           ),
                           showlegend=True,
                           name=clus + ' ' + name,
                           text=u.loc[selection_cond_bset].index
                           )
            )

            fig.add_trace(
                go.Scatter(mode='markers',
                           x=u.loc[selection_cond_aset, 'U1'],
                           y=u.loc[selection_cond_aset, 'U2'],
                           marker=dict(
                               color=u.loc[selection_cond_aset, 'cluster'].map(
                                   {'C1': palette[1], 'C2': palette[2], 'C3': palette[3], 'C4': palette[4], 'C5': palette[5]}
                               ),
                               size=10,
                               opacity=1.0,
                           ),
                           showlegend=True,
                           name=clus + ' ' + name2,
                           text=u.loc[selection_cond_aset].index
                           )
            )

        fig.update_layout(title_x=0.5, width=1000, height=800)
        fig.update_traces(marker={'size': 14})

        fig_json = fig.to_json()

        base_fn = 'umap_plotly_' + name + '_' + name2 + '.json'
        with open('../plots/' + base_fn, 'w+') as f:
            json.dump(fig_json, f, indent=4)

        with open('/Users/twood/Desktop/dlbcl-classifier-portal/files/plots/' + base_fn, 'w+') as f:
            json.dump(fig_json, f, indent=4)

for bset, name, bpreds in base_sets:

    bset_icm_json = make_icomut_json(bset.drop(['cluster', 'set', 'confidence']), qv, bset.loc[['cluster', 'set', 'confidence']])

    bpreds['set'] = name

    reduced_df = fd.construct_reduced_winning_version(bset.drop(['set', 'cluster', 'confidence'], axis=0), add_missing_features=True)
    pred_df = bpreds

    pred_df = pred_df.sort_values(by=['set', 'PredictedCluster', 'Confidence'], ascending=[True, True, True])

    n_comp = 2
    fit = umap.UMAP(n_components=n_comp, n_neighbors=15, min_dist=0.10, metric='euclidean', random_state=seed)
    u = pd.DataFrame(fit.fit_transform(reduced_df))
    u.index = reduced_df.index
    u.columns = ['U' + str(i) for i in range(1, n_comp + 1)]

    u['sample'] = u.index
    u['cluster'] = bset.loc['cluster']

    centroids = {g[0]: (g[1]['U1'].mean(), g[1]['U2'].mean()) for g in u.groupby('cluster')}
    u = get_umap_distances(u, centroids)

    pred_df.loc[:, 'U1'] = u.loc[pred_df.index, 'U1']
    pred_df.loc[:, 'U2'] = u.loc[pred_df.index, 'U2']
    pred_df.index.name = 'sample'

    pred_df.to_csv('/Users/twood/Desktop/dlbcl-classifier-portal/files/classifications/' + name + '_' + 'noset2' + '_predictions.tsv', sep='\t', float_format='%.4f')

    fn = '/Users/twood/Desktop/dlbcl-classifier-portal/files/data/dlbcl_' + name + '_' + 'noset2' + '_data.json'
    fn2 = '../portal_icomut_jsons/dlbcl_' + name + '_' + 'noset2' + '_data.json'
    with open(fn, 'w+') as f:
        json.dump(bset_icm_json, f, indent=4)
    with open(fn2, 'w+') as f:
        json.dump(bset_icm_json, f, indent=4)

    with open('/Users/twood/Desktop/dlbcl-classifier-portal/files/config/config.json', 'r') as f:
        cfg_file = json.load(f)

    del cfg_file['samplesToHighlight']
    del cfg_file['highlightColor']
    cfg_file['panels'][0]['label'] = "Cluster"

    set_panel = {
                    "type": "strip",
                    "id": "Set",
                    "label": "Set",
                    "dataPrefix": "set",
                    "row": 1,
                    "colors": [
                        "#000000",
                        "#f5f242"
                    ],
                    "plots": {
                        "center": {
                            "yAxisTitle": "Set"
                        }
                    }
                }
    confidence_panel = {
                            "type": "continuousstrip",
                            "id": "Confidence",
                            "label": "Confidence",
                            "dataPrefix": "confidence",
                            "row": 2,
                            "colors": [
                                "#ff1900",
                                "#00fa08"
                            ],
                            "plots": {
                                "center": {
                                    "yAxisTitle": "Confidence",
                                    "dataOptions": {"colorRange": {"min": 0.20, "max": 1.0}}
                                }
                            }
                        }

    cfg_file['panels'].insert(1, set_panel)
    cfg_file['panels'].insert(2, confidence_panel)
    cfg_file['panels'][3]['row'] = 3

    if name == 'Chapuy':
        cfg_file['panels'][0]['colors'] = ["#17becf", "#d62728", "#ff7f0e", "#2ca02c", "#9467bd"]

    if name == 'Schmitz':
        cfg_file['panels'][0]['colors'] = ["#9467bd", "#ff7f0e", "#d62728", "#2ca02c", "#17becf",]

    for i in range(0, len(cfg_file['panels'])):
        cfg_file['panels'][i]['displayLegend'] = True

    cfg_file['panels'][0]['displayLegend'] = False

    cfg_file['autoAddPadding'] = True

    cfg_file['initialSort'] = [{"panel": "Cluster", "order": 1},
                               {"panel": "Set", "order": 1},
                               {"panel": "Confidence", "order": 1}]

    fn = '/Users/twood/Desktop/dlbcl-classifier-portal/files/config/config_' + name + '_' + 'noset2' + '.json'
    fn2 = '../portal_config_jsons/config_' + name + '_' + 'noset2' + '.json'
    with open(fn, 'w+') as f:
        json.dump(cfg_file, f, indent=4)
    with open(fn2, 'w+') as f:
        json.dump(cfg_file, f, indent=4)


    #########
    # UMAPs #
    #########

    samples = bset.columns

    fig = go.Figure()

    for clus in ['C1', 'C2', 'C3', 'C4', 'C5']:
        selection_cond_bset = u.index.isin(samples) & (u['cluster'] == clus)

        fig.add_trace(
            go.Scatter(mode='markers',
                       x=u.loc[selection_cond_bset, 'U1'],
                       y=u.loc[selection_cond_bset, 'U2'],
                       marker=dict(
                           color=u.loc[selection_cond_bset, 'cluster'].map(
                               {'C1': palette[1], 'C2': palette[2], 'C3': palette[3], 'C4': palette[4], 'C5': palette[5]}
                           ),
                           size=10,
                           opacity=1.0,
                       ),
                       showlegend=True,
                       name=clus + ' ' + name,
                       text=u.loc[selection_cond_bset].index
                       )
        )

    fig.update_layout(title_x=0.5, width=1000, height=800)
    fig.update_traces(marker={'size': 14})

    fig_json = fig.to_json()

    base_fn = 'umap_plotly_' + name + '_' + 'noset2' + '.json'
    with open('../plots/' + base_fn, 'w+') as f:
        json.dump(fig_json, f, indent=4)

    with open('/Users/twood/Desktop/dlbcl-classifier-portal/files/plots/' + base_fn, 'w+') as f:
        json.dump(fig_json, f, indent=4)

