import pandas as pd
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import random

seed = 1234
random.seed(seed)
np.random.seed(seed)


umap = pd.read_csv('../data_tables/umap/umap_70conf.tsv', sep='\t', index_col=0)

preds = '../evaluation_validation_set/confidence_adjusted_tables/NN_reducedV3.2_removeN5_nfeatures21_pMax0.94248563.tsv'
preds_test = '../evaluation_test_set/NN_reducedV3.2_removeN5_nfeatures21_testsetEval.tsv'
gsm_file = '../data_tables/gsm/old_matrices/DLBCL_Staudt_Shipp_CL.for_classifier_training.classifier_subset.fix_sv.fix_ploidy.17-Aug-2022.txt'
training_set = list(pd.read_csv('../data_tables/train_test_sets/TrainingSet_550Subset_May2021.txt', sep='\t', header=None, index_col=0).index)
testing_set = list(pd.read_csv('../data_tables/train_test_sets/TestingSet_149Subset_May2021.txt', sep='\t', header=None, index_col=0).index)
sig_genes_file = '../data_tables/qval_dfs/fisher_exact_5x2_17-Aug-2022.combined.tsv'
coo_file = '../data_tables/phenotypes/'
set_file = '../data_tables/sample_sets/ShippStaudtSets.purity0.2.txt'

preds = pd.read_csv(preds, index_col=0, sep='\t')
preds_test = pd.read_csv(preds_test, index_col=0, sep='\t')
preds['Confidence'] = preds.max(axis=1)
confidences = pd.concat([preds_test['Confidence'], preds['Confidence']])
confidences = confidences.loc[umap.index]

sample_cohorts = pd.read_csv(set_file, sep='\t', index_col=0)

gsm = pd.read_csv(gsm_file, sep='\t', index_col=0).T
sig_genes = pd.read_csv(sig_genes_file, sep='\t', index_col=0)
sig_genes = sig_genes.loc[sig_genes['q'] <= 0.10].index

gsm_b = gsm[sig_genes].astype(float).astype(int) != 0
gsm_b['COO'] = gsm['COO']
gsm_b['PLOIDY'] = gsm['PLOIDY'].astype(float).round(2)

gsm_b = gsm_b.loc[umap.index]

full_df = pd.concat([umap, gsm_b], axis=1)

full_df['set'] = 'Train'
full_df.loc[full_df.index.isin(testing_set), 'set'] = 'Test'

full_df['cohort'] = sample_cohorts.loc[full_df.index, 'cohort']
full_df['Confidence'] = confidences


def annotate_umap_genes():
    for g in sig_genes:
        t_s = full_df.loc[full_df[g]]
        f_s = full_df.loc[~full_df[g]]

        plt.scatter(t_s['u1'], t_s['u2'], c='#d13328', label=g)
        plt.scatter(f_s['u1'], f_s['u2'], c='black', label='No ' + g)

        plt.legend()

        plt.savefig('../plots/umap/umap_annotations/genes/umap_'+g + '.pdf')
        plt.clf()


def annotate_umap_coo():
    abc_s = full_df.loc[full_df['COO'] == 'ABC']
    gcb_s = full_df.loc[full_df['COO'] == 'GCB']
    unc_s = full_df.loc[full_df['COO'] == 'UNC']

    plt.scatter(abc_s['u1'], abc_s['u2'], c='#d13328', label='ABC')
    plt.scatter(gcb_s['u1'], gcb_s['u2'], c='#194bff', label='GCB')
    plt.scatter(unc_s['u1'], unc_s['u2'], c='#f5fc72', label='Unclassified')

    plt.legend()
    plt.savefig('../plots/umap/umap_annotations/umap_coo.pdf')
    plt.clf()


def annotate_umap_traintest():
    train_s = full_df.loc[full_df['set'] == 'Train']
    test_s = full_df.loc[full_df['set'] == 'Test']

    plt.scatter(train_s['u1'], train_s['u2'], c='black', label='Train')
    plt.scatter(test_s['u1'], test_s['u2'], c='#d13328', label='Test')

    plt.legend()
    plt.savefig('../plots/umap/umap_annotations/umap_traintest.pdf')
    plt.clf()


def annotate_umap_cohort():
    shipp_s = full_df.loc[full_df['cohort'] == 'Shipp']
    staudt_s = full_df.loc[(full_df['cohort'] == 'Staudt') | (full_df['cohort'] == 'TCGA')]

    plt.scatter(shipp_s['u1'], shipp_s['u2'], c='#d13328', label='Chapuy et al.')
    plt.scatter(staudt_s['u1'], staudt_s['u2'], c='black', label='Schmitz et al.')

    plt.legend()
    plt.savefig('../plots/umap/umap_annotations/umap_cohort.pdf')
    plt.clf()


def annotate_umap_confidence():
    plt.scatter(full_df['u1'], full_df['u2'], c=full_df['Confidence'], cmap='Reds')
    plt.colorbar()
    plt.clim((0.70, 1.0))
    plt.savefig('../plots/umap/umap_annotations/umap_confidences.pdf')
    plt.clf()


annotate_umap_confidence()
