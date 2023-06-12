import pandas as pd
from matplotlib import pyplot as plt
import numpy as np
import random
import format_data as fd
from matplotlib.colors import LinearSegmentedColormap
import seaborn as sns

seed = 1234
random.seed(seed)
np.random.seed(seed)

colors = [(0, 0, 0), (1, 0, 0)]  # first color is black, last is red
cm = LinearSegmentedColormap.from_list("Custom", colors, N=20)


umap = pd.read_csv('../data_tables/umap/umap_70conf.tsv', sep='\t', index_col=0)
umap_full = pd.read_csv('../data_tables/umap/umap_all.tsv', sep='\t', index_col=0)
umap_full.columns = ['u1_full', 'u2_full']

preds = "../evaluation_validation_set/confidence_adjusted_tables/NN_reducedV3.4_removeN5_nfeatures21_pMax0.93856484.tsv"
preds_test = '../evaluation_test_set/NN_reducedV3.4_removeN5_nfeatures21_testsetEval.tsv'
gsm_file = '../data_tables/gsm/DLBCL.699.fullGSM.Sep_23_2022.tsv'
training_set = list(pd.read_csv('../data_tables/train_test_sets/TrainingSet_550Subset_May2021.txt', sep='\t', header=None, index_col=0).index)
testing_set = list(pd.read_csv('../data_tables/train_test_sets/TestingSet_149Subset_May2021.txt', sep='\t', header=None, index_col=0).index)
sig_genes_file = '../data_tables/qval_dfs/fisher_exact_5x2.Sep_23_2022.combined.tsv'
set_file = '../data_tables/sample_sets/ShippStaudtSets.purity0.2.txt'

preds = pd.read_csv(preds, index_col=0, sep='\t')
preds_test = pd.read_csv(preds_test, index_col=0, sep='\t')
preds['Confidence'] = preds.max(axis=1)
confidences2 = pd.concat([preds_test['Confidence'], preds['Confidence']])
confidences = confidences2.loc[umap.index]

sample_cohorts = pd.read_csv(set_file, sep='\t', index_col=0)

gsm = pd.read_csv(gsm_file, sep='\t', index_col=0).T
sig_genes = pd.read_csv(sig_genes_file, sep='\t', index_col=0)
sig_genes = sig_genes.loc[sig_genes['q'] <= 0.10].index

gsm_b = gsm[sig_genes].astype(float).astype(int) != 0
gsm_b['COO'] = gsm['COO']
gsm_b['PLOIDY'] = gsm['PLOIDY'].astype(float).round(2)

gsm_b = gsm_b.loc[umap.index]

full_df = pd.concat([umap, gsm_b], axis=1)
full_df2 = pd.concat([umap_full, gsm_b], axis=1)

full_df['set'] = 'Train'
full_df.loc[full_df.index.isin(testing_set), 'set'] = 'Test'
full_df2['set'] = 'Train'
full_df2.loc[full_df2.index.isin(testing_set), 'set'] = 'Test'

full_df['cohort'] = sample_cohorts.loc[full_df.index, 'cohort']
full_df['Confidence'] = confidences
full_df2['cohort'] = sample_cohorts.loc[full_df2.index, 'cohort']
full_df2['Confidence'] = confidences2

reduced_df2 = fd.construct_reduced_winning_version(gsm)
reduced_df = reduced_df2.loc[full_df.index]


def annotate_umap_genes():
    for g in sig_genes:
        plt.figure(figsize=(14, 12))
        t_s = full_df.loc[full_df[g]]
        f_s = full_df.loc[~full_df[g]]

        plt.scatter(t_s['u1'], t_s['u2'], c='#d13328', label=g, s=60)
        plt.scatter(f_s['u1'], f_s['u2'], c='black', label='No ' + g, s=60)

        plt.legend()

        plt.savefig('../plots/umap/umap_annotations/genes/umap_'+g + '.pdf')
        plt.clf()
        plt.close()


def annotate_umap_coo():
    plt.figure(figsize=(14, 12))
    abc_s = full_df.loc[full_df['COO'] == 'ABC']
    gcb_s = full_df.loc[full_df['COO'] == 'GCB']
    unc_s = full_df.loc[full_df['COO'] == 'UNC']

    plt.scatter(abc_s['u1'], abc_s['u2'], c='#d13328', label='ABC', s=60)
    plt.scatter(gcb_s['u1'], gcb_s['u2'], c='#194bff', label='GCB', s=60)
    plt.scatter(unc_s['u1'], unc_s['u2'], c='#f5fc72', label='Unclassified', s=60)

    plt.legend()
    plt.savefig('../plots/umap/umap_annotations/umap_coo.pdf')
    plt.clf()


def annotate_umap_traintest():
    plt.figure(figsize=(14, 12))
    train_s = full_df.loc[full_df['set'] == 'Train']
    test_s = full_df.loc[full_df['set'] == 'Test']

    plt.scatter(train_s['u1'], train_s['u2'], c='black', label='Train', s=60)
    plt.scatter(test_s['u1'], test_s['u2'], c='#d13328', label='Test', s=60)

    plt.legend()
    plt.savefig('../plots/umap/umap_annotations/umap_traintest.pdf')
    plt.clf()


def annotate_umap_cohort():
    plt.figure(figsize=(14, 12))
    shipp_s = full_df.loc[full_df['cohort'] == 'Shipp']
    staudt_s = full_df.loc[(full_df['cohort'] == 'Staudt') | (full_df['cohort'] == 'TCGA')]

    plt.scatter(shipp_s['u1'], shipp_s['u2'], c='#d13328', label='Chapuy et al.', s=60)
    plt.scatter(staudt_s['u1'], staudt_s['u2'], c='black', label='Schmitz et al.', s=60)

    plt.legend()
    plt.savefig('../plots/umap/umap_annotations/umap_cohort.pdf')
    plt.clf()


def annotate_umap_confidence():
    colors = [(1, 0, 0), (0, 0, 0)]
    cm = LinearSegmentedColormap.from_list("Custom", colors, N=20)
    plt.figure(figsize=(14, 12))
    plt.scatter(full_df['u1'], full_df['u2'], c=full_df['Confidence'], cmap=cm, s=60)
    cbar = plt.colorbar()
    cbar.set_label('Confidence', rotation=270, size=50, labelpad=70)
    plt.clim((0.70, 1.0))
    plt.savefig('../plots/umap/umap_annotations/umap_confidences.pdf')
    plt.clf()

    plt.figure(figsize=(14, 12))
    plt.scatter(full_df2['u1_full'], full_df2['u2_full'], c=full_df2['Confidence'], cmap=cm, s=60)
    cbar = plt.colorbar()
    cbar.set_label('Confidence', rotation=270, size=50, labelpad=70)
    plt.clim((0.2, 1.0))
    plt.savefig('../plots/umap/umap_annotations/umap_confidences_full.pdf')
    plt.clf()


def umap_reduced_features():
    c1_features = ['BCL6_ALT', 'NOTCH2_vec', 'M88O_vec', 'C1_vec4', 'CD70_vec']
    c2_features = ['TP53_biallelic', 'X21Q_AMP', 'Sum_C2_ARM', 'Sum_C2_FOCAL']
    c3_features = ['BCL2_combined', 'CREBBP_vec', 'GNA13_vec', 'PTEN', 'SV_MYC']
    c4_features = ['Hist_comp', 'SGK1_vec', 'DUSP2_vec']
    c5_features = ['TBL1XR1_vec', 'MYD88_L265P_CD79B', 'Sum_C5_CNA']
    misc1 = ['CN_2P16_1_AMP']

    sum_c1 = reduced_df[c1_features].sum(axis=1)
    sum_c2 = reduced_df[c2_features].sum(axis=1)
    sum_c3 = reduced_df[c3_features].sum(axis=1)
    sum_c4 = reduced_df[c4_features].sum(axis=1)
    sum_c5 = reduced_df[c5_features].sum(axis=1)
    sum_misc1 = reduced_df[misc1].sum(axis=1)

    plot_list = [(sum_c1, 'C1'), (sum_c2, 'C2'), (sum_c3, 'C3'),
                 (sum_c4, 'C4'), (sum_c5, 'C5'), (sum_misc1, '2P16.1.AMP')]

    pallete = {1: sns.color_palette()[4],  # purple, C1
              2: sns.color_palette()[9],  # blue, C2
              3: sns.color_palette()[1],  # orange, C3
              4: sns.color_palette()[2],  # green, C4
              5: sns.color_palette()[3],  # red, C5
              6: (0, 0, 1)
              }

    for i, p in enumerate(plot_list):
        clus_color = [(1, 1, 1), pallete[i + 1]]
        plt.figure(figsize=(14, 12))
        cm_clus = LinearSegmentedColormap.from_list("Custom", clus_color, N=20)
        plt.scatter(full_df['u1'], full_df['u2'], c=p[0], cmap=cm_clus, s=60)
        cbar = plt.colorbar()
        cbar.set_label('Reduced Cluster Sum', rotation=270, size=50, labelpad=70)
        plt.title(p[1], size=50)
        plt.savefig('../plots/umap/umap_annotations/umap_reduced_' + p[1] + '.pdf')
        plt.clf()
        plt.close()


#annotate_umap_genes()
#annotate_umap_coo()
#annotate_umap_traintest()
#annotate_umap_cohort()
annotate_umap_confidence()
#umap_reduced_features()
