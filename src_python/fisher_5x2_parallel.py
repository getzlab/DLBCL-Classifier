import os
if "R_HOME" not in os.environ:
    os.environ['R_HOME'] = '/Library/Frameworks/R.framework/Resources/'
import numpy as np
import rpy2.robjects.numpy2ri
import rpy2.robjects as R
from rpy2.robjects.packages import importr
rpy2.robjects.numpy2ri.activate()
import time
import multiprocessing as mp
import pandas as pd
import statsmodels.stats.multitest as sm

R.r('set.seed')(1)

R_STATS = importr('stats')
BLACKLIST = ['MUC6', 'HIST1H2BK', 'HIST2H2BE', 'OR51B6', 'OR10V1']
VERSION_DATE = 'Sep_23_2022'

MODE = ['combined', 'staudtonly', 'shipponly']
MODE = MODE[0]


def fisher_exact_2x2(matrix, alt):
    return R_STATS.fisher_test(matrix, alternative=alt)


def fisher_exact_5x2(matrix, numiter=100000):
    return R_STATS.fisher_test(matrix, simulate_p_value=True, B=numiter)


labels = '../data_tables/clustering_labels/GSM699_cluster_Sep_23_2022.bestclus.remapped.txt'
training_file = '../data_tables/train_test_sets/TrainingSet_550Subset_May2021.txt'
gsm_file = '../data_tables/gsm/DLBCL.699.fullGSM.Sep_23_2022.tsv'

cohorts = pd.read_csv('../data_tables/sample_sets/ShippStaudtSets.purity0.2.txt', index_col=0, sep='\t')

clusters = pd.read_csv(labels, sep='\t', index_col=0)
gsm = pd.read_csv(gsm_file, sep='\t', index_col=0)
training_samples = list(pd.read_csv(training_file, sep='\t', index_col=0, header=None).index)
staudt_samples = cohorts.loc[cohorts['cohort'] != 'Shipp'].index
shipp_samples = cohorts.loc[cohorts['cohort'] == 'Shipp'].index

if MODE == 'staudtonly':
    gsm = gsm[gsm.columns.intersection(staudt_samples)]
    clusters = pd.read_csv('../data_tables/clustering_labels/GSM699_cluster_Sep_23_2022_STAUDT.bestclus.txt', sep='\t', index_col=0)
if MODE == 'shipponly':
    gsm = gsm[gsm.columns.intersection(shipp_samples)]
    clusters = pd.read_csv('../data_tables/clustering_labels/GSM699_cluster_Sep_23_2022_SHIPP.bestclus.remapped.txt', sep='\t', index_col=0)

gsm.loc['PLOIDY'] = (gsm.loc['PLOIDY'].astype(float) > 2.5).astype(int)
gsm.loc['COO_ABC'] = gsm.loc['COO'].map({'ABC': 1, 'GCB': 0, 'UNC': 0, 'na': 'na'})
gsm.loc['COO_GCB'] = gsm.loc['COO'].map({'ABC': 0, 'GCB': 1, 'UNC': 0, 'na': 'na'})
gsm.loc['COO_UNC'] = gsm.loc['COO'].map({'ABC': 0, 'GCB': 0, 'UNC': 1, 'na': 'na'})
gsm = gsm.drop('COO')
gsm = gsm.drop('PURITY')
gsm = gsm.loc[~gsm.index.str.contains('CCF')]

rows = [i for i in gsm.index if i not in ['COO_ABC', 'COO_GCB', 'COO_UNC']]
for idx in rows:
    gsm.loc[idx] = gsm.loc[idx].astype(float).astype(int)

tables = []

gsm_train = gsm[gsm.columns.intersection(training_samples)]

# Minimize subset calls within the loop by pre-assigning variables

found_clusters = sorted(list(set(clusters['cluster'])))
sample_sets = []
cluster_dfs = []
for clus in found_clusters:
    ci_samples = list(clusters.loc[clusters['cluster'] == clus].index)
    ci_samples = [x for x in ci_samples if x in training_samples]

    ci_df = gsm_train[ci_samples]

    cluster_dfs.append(ci_df)
    sample_sets.append(ci_samples)

# Do 2x2 for COO and PLOIDY first
coo_ploidy = ['COO_ABC', 'COO_GCB', 'COO_UNC', 'PLOIDY']

cols = ['C' + str(i) + '_p_greater' for i in found_clusters] + \
       ['C' + str(i) + '_p_less' for i in found_clusters]
for i in found_clusters:
    cols += ['C' + str(i) + '_mut', 'C' + str(i) + '_wt', '!C' + str(i) + '_mut', '!C' + str(i) + '_wt']

stats_table_coo = np.array([[0] * len(cols)] * 4)
stats_table_coo = pd.DataFrame(stats_table_coo)
stats_table_coo.index = coo_ploidy

stats_table_coo.columns = cols

for event in coo_ploidy:
    print(event)
    for i in range(len(sample_sets)):
        samples = sample_sets[i]
        cluster = str(i+1)
        clus_df = cluster_dfs[i]
        non_clus_df = gsm_train.drop(samples, axis=1)

        clus_mt = (clus_df.loc[event] == 1).sum()
        clus_wt = (clus_df.loc[event] == 0).sum()
        non_mt = (non_clus_df.loc[event] == 1).sum()
        non_wt = (non_clus_df.loc[event] == 0).sum()

        m = np.array([[0, 0]] * 2)
        m[0][0] = clus_mt
        m[0][1] = clus_wt
        m[1][0] = non_mt
        m[1][1] = non_wt

        p_g = fisher_exact_2x2(m, 'greater')[0][0]
        p_l = fisher_exact_2x2(m, 'less')[0][0]

        stats_table_coo.loc[event, 'C' + cluster + '_p_greater'] = p_g
        stats_table_coo.loc[event, 'C' + cluster + '_p_less'] = p_l
        stats_table_coo.loc[event, 'C' + cluster + '_mut'] = clus_mt
        stats_table_coo.loc[event, 'C' + cluster + '_wt'] = clus_wt
        stats_table_coo.loc[event, '!C' + cluster + '_mut'] = non_mt
        stats_table_coo.loc[event, '!C' + cluster + '_wt'] = non_wt


stats_table_coo.to_csv('../data_tables/qval_dfs/fisher_exact_2x2.coo_ploidy.' + VERSION_DATE + '.' + MODE + '.tsv', sep='\t')

# 5x2 tests

cols = ['p'] + ['C' + str(i) + '_nf' for i in found_clusters] + ['C' + str(i) + '_f' for i in found_clusters]

for i in found_clusters:
    cols += ['C' + str(i) + '_mut', 'C' + str(i) + '_wt']

cols += ['C' + str(i) + '_sum' for i in found_clusters]

stats_table = np.array([[0] * len(cols)] * len(gsm_train.index))
stats_table = pd.DataFrame(stats_table)
stats_table.index = gsm_train.index
stats_table.columns = cols
stats_table['p'] = -1

for gene in gsm_train.index:
    print(gene)
    m = np.array([[0, 0]] * 5)
    for i in range(len(cluster_dfs)):
        ci_df = cluster_dfs[i]
        clus = 'C' + str(i + 1)
        if 'COO' in gene:
            n_mut = (ci_df.loc[gene] == 1).sum()
            n_wt = (ci_df.loc[gene] == 0).sum()
        else:
            n_mut = (ci_df.loc[gene] > 0).sum()
            n_wt = (ci_df.loc[gene] == 0).sum()
            stats_table.loc[gene, clus + '_sum'] = ci_df.loc[gene].sum()
        m[i][0] = n_mut
        m[i][1] = n_wt
        stats_table.loc[gene, clus+'_mut'] = n_mut
        stats_table.loc[gene, clus+'_wt'] = n_wt

    if sum(m[:, 0]) != 0:
        p = fisher_exact_5x2(m)[0][0]
        stats_table.loc[gene, 'p'] = p

denom = len(gsm_train.columns)
print('Denom =', denom)

clus_names = ['C' + str(i) for i in found_clusters]
for f_c in clus_names:
    stats_table[f_c + '_f'] = stats_table[f_c + '_mut'] / denom

denom = 0
for f_c in clus_names:
    denom += stats_table[f_c + '_f']

for f_c in clus_names:
    stats_table[f_c + '_nf'] = stats_table[f_c + '_f'] / denom

# Drop genes in the BLACKLIST or genes < 1% frequency.
# Drop COO rows and PLOIDY as well, since we calculated these separately

o_f = 0
for f_c in clus_names:
    o_f += stats_table[f_c + '_f']

stats_table.insert(1, 'overall_frequency', o_f)

stats_table = stats_table.drop(coo_ploidy)
stats_table = stats_table.loc[stats_table['overall_frequency'] >= 0.01]

# q value calculation

stats_table.insert(0, 'q', sm.multipletests(stats_table['p'], method='fdr_bh')[1])

nf_cols = ['C' + str(i) + '_nf' for i in found_clusters]
mapping = {'C' + str(i) + '_nf': 'C' + str(i) for i in found_clusters}
stats_table.insert(0, 'cluster', stats_table[nf_cols].idxmax(axis=1).map(mapping))

clus_order = []
for c in clus_names:
    c_o = list(stats_table.loc[stats_table['cluster'] == c].sort_values(by=['q', c + '_f'], ascending=[True, False]).index)
    clus_order += c_o

stats_table = stats_table.loc[clus_order]
stats_table.to_csv('../data_tables/qval_dfs/fisher_exact_5x2.' + VERSION_DATE + '.' + MODE + '.tsv', sep='\t')


# How to parallelize for the future
# Step 1: Init multiprocessing.Pool()
# pool = mp.Pool(mp.cpu_count())
#
# # Step 2: `pool.apply` the `howmany_within_range()`
# results = [pool.apply(howmany_within_range, args=(row, 4, 8)) for row in data]
#
# # Step 3: Don't forget to close
# pool.close()
#
# print(results[:10])

# pool = mp.Pool(mp.cpu_count())
#
#
#
#
# r_stats = importr('stats')
# m = np.array([[4, 100], [500, 2], [30, 30], [1, 1], [10, 30]])
# start = time.time()
# res = r_stats.fisher_test(m, simulate_p_value=True, B=1000000)
# end = time.time()
#
# print(res)