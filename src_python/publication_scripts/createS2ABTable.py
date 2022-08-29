import os
if "R_HOME" not in os.environ:
    os.environ['R_HOME'] = '/Library/Frameworks/R.framework/Resources/'
import pandas as pd
import rpy2.robjects.numpy2ri
import rpy2.robjects as R
from rpy2.robjects.packages import importr
import scipy.stats as ss
import numpy as np
from matplotlib import pyplot as plt
from lifelines import KaplanMeierFitter
from lifelines.utils import median_survival_times
from lifelines.statistics import pairwise_logrank_test
from decimal import Decimal

rpy2.robjects.numpy2ri.activate()

R.r('set.seed')(1)
R_STATS = importr('stats')


def pg(vec, x):
    return str(np.round((vec[x] / sum(vec)) * 100, 1))


def fisher_exact_2x2(matrix, alt='two.sided'):
    return R_STATS.fisher_test(matrix, alternative=alt)


def fisher_exact_mxn(matrix, numiter=1000000):
    return R_STATS.fisher_test(matrix, simulate_p_value=True, B=numiter)


s1_t = pd.read_csv('../../data_tables/tableS1_classifier_merged.tsv', sep='\t', index_col=0)
print(s1_t.columns)
# Summary table for S2A
# 		         Full cohort	Schmitz et al		Chapuy et al 		p value         test used
# Sample (N)
# F/M (N) - fisher exact
# ABC/GCB/unclass - 3x2 fisher
# TO vs not - fisher
# r-chop - fisher
# Age - mann-whitney
# IPI - fisher 5x2 test
# PFS years (convert to months), censoring is PFS.status - reversed log rank
# OFS (os.time), censoring is OS.status - reversed log rank

table_s2a = pd.DataFrame(columns=['Full Cohort', 'Schmitz et al.', 'Chapuy et al.', 'p value', 'test'])
full_counts = s1_t['Cohort'].value_counts()
table_s2a.loc['Samples (N)'] = ['699',
                                full_counts['Schmitz et al.'],
                                full_counts['Chapuy et al.'],
                                '', '']

# fisher test for Sex (3x2 fisher)
sex_counts = s1_t['Gender'].value_counts()
sex_f = ';'.join([x + '=' + str(sex_counts[x]) for x in sex_counts.index])
sex_a = s1_t.groupby('Cohort')['Gender'].value_counts()
sex_ch = sex_a['Chapuy et al.']
sex_sc = sex_a['Schmitz et al.']

sex_t = np.matrix([sex_ch.values, sex_sc.values])
sex_p = fisher_exact_2x2(sex_t)[0][0]


table_s2a.loc['Sex'] = [sex_f,
                        ';'.join([x + '=' + str(sex_sc[x]) + '(' + pg(sex_sc, x) + '%)' for x in sex_sc.index]),
                        ';'.join([x + '=' + str(sex_ch[x]) + '(' + pg(sex_ch, x) + '%)' for x in sex_ch.index]),
                        sex_p, 'Fisher Exact (2x2)']

# fisher test for COO (3x2 fisher)
coo_counts = s1_t['COO'].value_counts()
coo_f = ';'.join([x + '=' + str(coo_counts[x]) for x in coo_counts.index])
coo_a = s1_t.groupby('Cohort')['COO'].value_counts()
coo_ch = coo_a['Chapuy et al.']
coo_sc = coo_a['Schmitz et al.']

coo_t = np.matrix([coo_ch.values, coo_sc.values])
coo_p = fisher_exact_mxn(coo_t)[0][0]

table_s2a.loc['COO'] = [coo_f,
                        ';'.join([x + '=' + str(coo_sc[x]) + '(' + pg(coo_sc, x) + '%)' for x in coo_sc.index]),
                        ';'.join([x + '=' + str(coo_ch[x]) + '(' + pg(coo_ch, x) + '%)' for x in coo_ch.index]),
                        coo_p, 'Fisher Exact (3x2)']

# fisher test for TO vs TN

sp_counts = s1_t['Pair Status (Tumor Only/Normal)'].value_counts()
sp_f = ';'.join([x + '=' + str(sp_counts[x]) for x in sp_counts.index])
sp_a = s1_t.groupby('Cohort')['Pair Status (Tumor Only/Normal)'].value_counts()
sp_ch = sp_a['Chapuy et al.']
sp_sc = sp_a['Schmitz et al.']

sp_t = np.matrix([sp_ch.values, sp_sc.values])
sp_p = fisher_exact_2x2(sp_t)[0][0]

table_s2a.loc['Pair Status (TN/TO)'] = [sp_f,
                                        ';'.join([x + '=' + str(sp_sc[x]) + '(' + pg(sp_sc, x) + '%)' for x in sp_sc.index]),
                                        ';'.join([x + '=' + str(sp_ch[x]) + '(' + pg(sp_ch, x) + '%)' for x in sp_ch.index]),
                                        sp_p, 'Fisher Exact (2x2)']

# fisher test for r-chop
rc_counts = s1_t['R-CHOp-like\nChemo'].value_counts()
rc_f = ';'.join([x + '=' + str(rc_counts[x]) for x in rc_counts.index])
rc_a = s1_t.groupby('Cohort')['R-CHOp-like\nChemo'].value_counts()
rc_ch = rc_a['Chapuy et al.']
rc_sc = rc_a['Schmitz et al.']

rc_t = np.matrix([rc_ch.values, rc_sc.values])
rc_p = fisher_exact_2x2(rc_t)[0][0]

table_s2a.loc['R-CHOp-like Chemo'] = [rc_f,
                                      ';'.join([x + '=' + str(rc_sc[x]) + '(' + pg(rc_sc, x) + '%)' for x in rc_sc.index]),
                                      ';'.join([x + '=' + str(rc_ch[x]) + '(' + pg(rc_ch, x) + '%)' for x in rc_ch.index]),
                                      rc_p, 'Fisher Exact (2x2)']



# Mann-Whitney U test for age
age_med = s1_t['Age-at first diagnosis'].median()
age_a_med = s1_t.groupby('Cohort')['Age-at first diagnosis'].median()
age_ch_med = age_a_med['Chapuy et al.']
age_sc_med = age_a_med['Schmitz et al.']

ch_ages = s1_t.groupby('Cohort')['Age-at first diagnosis'].get_group('Chapuy et al.').dropna()
sc_ages = s1_t.groupby('Cohort')['Age-at first diagnosis'].get_group('Schmitz et al.').dropna()

sc_ages = sc_ages[sc_ages != 0.0]

_, age_p = ss.mannwhitneyu(ch_ages, sc_ages, method="asymptotic")

table_s2a.loc['Age-at first diagnosis (Median)'] = [age_med,
                                                    age_sc_med,
                                                    age_ch_med,
                                                    age_p, 'Mann-Whitney U']


# IPI - fisher 6x2 test
ipi_counts = s1_t['IPI'].value_counts()
ipi_f = ';'.join([str(int(x)) + '=' + str(ipi_counts[x]) + '(' + pg(ipi_counts, x) + '%)' for x in sorted(ipi_counts.index)])
ipi_a = s1_t.groupby('Cohort')['IPI'].value_counts()
ipi_ch = ipi_a['Chapuy et al.'].sort_index()
ipi_sc = ipi_a['Schmitz et al.'].sort_index()

ipi_t = np.matrix([ipi_ch, ipi_sc])
ipi_p = fisher_exact_mxn(ipi_t)[0][0]

table_s2a.loc['IPI'] = [ipi_f,
                        ';'.join([str(int(x)) + '=' + str(ipi_sc[x]) + '(' + pg(ipi_sc, x) + '%)' for x in sorted(ipi_sc.index)]),
                        ';'.join([str(int(x)) + '=' + str(ipi_ch[x]) + '(' + pg(ipi_ch, x) + '%)' for x in sorted(ipi_ch.index)]),
                        ipi_p, 'Fisher Exact (6x2)']

# PFS years

tmp = s1_t.loc[~s1_t['PFS-years'].isna()]
log_rank = pairwise_logrank_test(tmp['PFS-years'], tmp['PFS_STAT'], tmp['Cohort'])

pfs_med = np.round(s1_t['PFS-years'].median(), 2)
pfs_ch_med = np.round(s1_t.loc[s1_t['Cohort'] == 'Chapuy et al.', 'PFS-years'].median(), 2)
pfs_sc_med = np.round(s1_t.loc[s1_t['Cohort'] == 'Schmitz et al.', 'PFS-years'].median(), 2)

table_s2a.loc['PFS (years, Median)'] = [pfs_med, pfs_sc_med, pfs_ch_med, log_rank.summary['p'].values[0], 'Log Rank']

# OS
tmp = s1_t.loc[~s1_t['OS.time'].isna()]
log_rank = pairwise_logrank_test(tmp['OS.time'], tmp['OS.status (1=dead)'], tmp['Cohort'])

os_med = np.round(s1_t['OS.time'].median(), 2)
os_ch_med = np.round(s1_t.loc[s1_t['Cohort'] == 'Chapuy et al.', 'OS.time'].median(), 2)
os_sc_med = np.round(s1_t.loc[s1_t['Cohort'] == 'Schmitz et al.', 'OS.time'].median(), 2)

table_s2a.loc['OS (years, Median)'] = [os_med, os_sc_med, os_ch_med, log_rank.summary['p'].values[0], 'Log Rank']

table_s2a['p value'] = table_s2a['p value'].apply(lambda x: x if x == '' else '%.2E' % Decimal(np.format_float_positional(x, precision=2, fractional=False)))
table_s2a.to_csv('../../data_tables/tableS2A_clinical.tsv', sep='\t')

# S2B
# Genetic Features:
# TO muts vs TN muts - mann-whitney
# Coverage - mann-whitney
# Purity - mann-whitney
# Ploidy - mann-whitney
# Mut Density - get from patientstats file from mutsig - mann-whitney

# Pair Status (Tumor Only/Normal)

table_s2b = pd.DataFrame(columns=['Full Cohort',
                                  'Schmitz et al.', 'Chapuy et al.', 'p value (Schmitz vs Chapuy)',
                                  'Tumor Only (TO)', 'Tumor Normal (TN)', 'p value (TO vs TN)', 'test'])

mut_a = s1_t.groupby('Cohort')['Number of Mutations']
mut_ch = mut_a.get_group('Chapuy et al.').median()
mut_sc = mut_a.get_group('Schmitz et al.').median()
mut_med = s1_t['Number of Mutations'].median()

ch_muts = s1_t.groupby('Cohort')['Number of Mutations'].get_group('Chapuy et al.').dropna()
sc_muts = s1_t.groupby('Cohort')['Number of Mutations'].get_group('Schmitz et al.').dropna()

mut_a = s1_t.groupby('Pair Status (Tumor Only/Normal)')['Number of Mutations']
mut_tn = mut_a.get_group('TN').median()
mut_to = mut_a.get_group('TO').median()

tn_muts = s1_t.groupby('Pair Status (Tumor Only/Normal)')['Number of Mutations'].get_group('TN').dropna()
to_muts = s1_t.groupby('Pair Status (Tumor Only/Normal)')['Number of Mutations'].get_group('TO').dropna()

_, muts_p = ss.mannwhitneyu(ch_muts, sc_muts, method="asymptotic")
_, muts_p2 = ss.mannwhitneyu(tn_muts, to_muts, method="asymptotic")

table_s2b.loc['Number of Mutations (Median)'] = [mut_med,
                                                 mut_sc, mut_ch, muts_p,
                                                 mut_to, mut_tn, muts_p2,
                                                 'Mann-Whitney U']

# Fraction deleted
fd_a = s1_t.groupby('Cohort')['Fraction Genome Deleted']
fd_ch = fd_a.get_group('Chapuy et al.').median()
fd_sc = fd_a.get_group('Schmitz et al.').median()
fd_med = s1_t['Fraction Genome Deleted'].median()

ch_fd = s1_t.groupby('Cohort')['Fraction Genome Deleted'].get_group('Chapuy et al.').dropna()
sc_fd = s1_t.groupby('Cohort')['Fraction Genome Deleted'].get_group('Schmitz et al.').dropna()

fd_a = s1_t.groupby('Pair Status (Tumor Only/Normal)')['Fraction Genome Deleted']
fd_tn = fd_a.get_group('TN').median()
fd_to = fd_a.get_group('TO').median()

tn_fd = s1_t.groupby('Pair Status (Tumor Only/Normal)')['Fraction Genome Deleted'].get_group('TN').dropna()
to_fd = s1_t.groupby('Pair Status (Tumor Only/Normal)')['Fraction Genome Deleted'].get_group('TO').dropna()

_, fd_p = ss.mannwhitneyu(ch_fd, sc_fd, method="asymptotic")
_, fd_p2 = ss.mannwhitneyu(tn_fd, to_fd, method="asymptotic")

table_s2b.loc['Fraction Genome Deleted (Median)'] = [fd_med,
                                                 fd_sc, fd_ch, fd_p,
                                                 fd_to, fd_tn, fd_p2,
                                                 'Mann-Whitney U']

# Fraction amplified
fa_a = s1_t.groupby('Cohort')['Fraction Genome Amplified']
fa_ch = fa_a.get_group('Chapuy et al.').median()
fa_sc = fa_a.get_group('Schmitz et al.').median()
fa_med = s1_t['Fraction Genome Amplified'].median()

ch_fa = s1_t.groupby('Cohort')['Fraction Genome Amplified'].get_group('Chapuy et al.').dropna()
sc_fa = s1_t.groupby('Cohort')['Fraction Genome Amplified'].get_group('Schmitz et al.').dropna()

fa_a = s1_t.groupby('Pair Status (Tumor Only/Normal)')['Fraction Genome Amplified']
fa_tn = fa_a.get_group('TN').median()
fa_to = fa_a.get_group('TO').median()

tn_fa = s1_t.groupby('Pair Status (Tumor Only/Normal)')['Fraction Genome Amplified'].get_group('TN').dropna()
to_fa = s1_t.groupby('Pair Status (Tumor Only/Normal)')['Fraction Genome Amplified'].get_group('TO').dropna()

_, fa_p = ss.mannwhitneyu(ch_fa, sc_fa, method="asymptotic")
_, fa_p2 = ss.mannwhitneyu(tn_fa, to_fa, method="asymptotic")

table_s2b.loc['Fraction Genome Amplified (Median)'] = [fa_med,
                                                 fa_sc, fa_ch, fa_p,
                                                 fa_to, fa_tn, fa_p2,
                                                 'Mann-Whitney U']

# Coverage

cov_a = s1_t.groupby('Cohort')['Mean Target Coverage']
cov_ch = np.round(cov_a.get_group('Chapuy et al.').median(), 2)
cov_sc = np.round(cov_a.get_group('Schmitz et al.').median(), 2)
cov_med = np.round(s1_t['Mean Target Coverage'].median(), 2)

ch_cov = s1_t.groupby('Cohort')['Mean Target Coverage'].get_group('Chapuy et al.').dropna()
sc_cov = s1_t.groupby('Cohort')['Mean Target Coverage'].get_group('Schmitz et al.').dropna()

cov_a = s1_t.groupby('Pair Status (Tumor Only/Normal)')['Mean Target Coverage']
cov_tn = np.round(cov_a.get_group('TN').median(), 2)
cov_to = np.round(cov_a.get_group('TO').median(), 2)

tn_cov = s1_t.groupby('Pair Status (Tumor Only/Normal)')['Mean Target Coverage'].get_group('TN').dropna()
to_cov = s1_t.groupby('Pair Status (Tumor Only/Normal)')['Mean Target Coverage'].get_group('TO').dropna()

_, cov_p = ss.mannwhitneyu(ch_cov, sc_cov, method="asymptotic")
_, cov_p2 = ss.mannwhitneyu(tn_cov, to_cov, method="asymptotic")

table_s2b.loc['Mean Target Coverage (Median)'] = [cov_med,
                                                  cov_sc, cov_ch, cov_p,
                                                  cov_to, cov_tn, cov_p2,
                                                  'Mann-Whitney U']


# Purity

pur_a = s1_t.groupby('Cohort')['Purity (ABSOLUTE)']
pur_ch = np.round(pur_a.get_group('Chapuy et al.').median(), 2)
pur_sc = np.round(pur_a.get_group('Schmitz et al.').median(), 2)
pur_med = np.round(s1_t['Purity (ABSOLUTE)'].median(), 2)

ch_pur = s1_t.groupby('Cohort')['Purity (ABSOLUTE)'].get_group('Chapuy et al.').dropna()
sc_pur = s1_t.groupby('Cohort')['Purity (ABSOLUTE)'].get_group('Schmitz et al.').dropna()

pur_a = s1_t.groupby('Pair Status (Tumor Only/Normal)')['Purity (ABSOLUTE)']
pur_tn = np.round(pur_a.get_group('TN').median(), 2)
pur_to = np.round(pur_a.get_group('TO').median(), 2)

tn_pur = s1_t.groupby('Pair Status (Tumor Only/Normal)')['Purity (ABSOLUTE)'].get_group('TN').dropna()
to_pur = s1_t.groupby('Pair Status (Tumor Only/Normal)')['Purity (ABSOLUTE)'].get_group('TO').dropna()

_, pur_p = ss.mannwhitneyu(ch_pur, sc_pur, method="asymptotic")
_, pur_p2 = ss.mannwhitneyu(tn_pur, to_pur, method="asymptotic")

table_s2b.loc['Purity (Median)'] = [pur_med,
                                    pur_sc,pur_ch, pur_p,
                                    pur_to, pur_tn, pur_p2,
                                    'Mann-Whitney U']

# ploidy

plo_a = s1_t.groupby('Cohort')['Ploidy (ABSOLUTE)']
plo_ch = np.round(plo_a.get_group('Chapuy et al.').median(), 2)
plo_sc = np.round(plo_a.get_group('Schmitz et al.').median(), 2)
plo_med = np.round(s1_t['Ploidy (ABSOLUTE)'].median(), 2)

ch_plo = s1_t.groupby('Cohort')['Ploidy (ABSOLUTE)'].get_group('Chapuy et al.').dropna()
sc_plo = s1_t.groupby('Cohort')['Ploidy (ABSOLUTE)'].get_group('Schmitz et al.').dropna()

plo_a = s1_t.groupby('Pair Status (Tumor Only/Normal)')['Ploidy (ABSOLUTE)']
plo_tn = np.round(plo_a.get_group('TN').median(), 2)
plo_to = np.round(plo_a.get_group('TO').median(), 2)

tn_plo = s1_t.groupby('Pair Status (Tumor Only/Normal)')['Ploidy (ABSOLUTE)'].get_group('TN').dropna()
to_plo = s1_t.groupby('Pair Status (Tumor Only/Normal)')['Ploidy (ABSOLUTE)'].get_group('TO').dropna()

_, plo_p = ss.mannwhitneyu(ch_plo, sc_plo, method="asymptotic")
_, plo_p2 = ss.mannwhitneyu(tn_plo, to_plo, method="asymptotic")

table_s2b.loc['Ploidy (Median)'] = [plo_med,
                                    plo_sc, plo_ch, plo_p,
                                    plo_to, plo_tn, plo_p2,
                                    'Mann-Whitney U']

# ploidy (boolean)

gd_a = s1_t.groupby('Cohort')['Ploidy (ABSOLUTE)']
gd_ch = (gd_a.get_group('Chapuy et al.') > 3.0).value_counts()
gd_sc = (gd_a.get_group('Schmitz et al.') > 3.0).value_counts()
gd_f = (s1_t['Ploidy (ABSOLUTE)'] > 3.0).value_counts()
gd_f = ';'.join([str(x) + '=' + str(gd_f[x]) for x in gd_f.index])

gd_a = s1_t.groupby('Pair Status (Tumor Only/Normal)')['Ploidy (ABSOLUTE)']
gd_to = (gd_a.get_group('TO') > 3.0).value_counts()
gd_tn = (gd_a.get_group('TN') > 3.0).value_counts()

gd_t = np.matrix([gd_ch.values, gd_sc.values])
gd_p = fisher_exact_2x2(gd_t)[0][0]

gd_t2 = np.matrix([gd_to.values, gd_tn.values])
gd_p2 = fisher_exact_2x2(gd_t2)[0][0]

table_s2b.loc['Genome Doubling (Fraction Doubled > 3.0)'] = [gd_f,
                                                             ';'.join([str(x) + '=' + str(gd_sc[x]) + '(' + pg(gd_sc, x) + '%)' for x in gd_sc.index]),
                                                             ';'.join([str(x) + '=' + str(gd_ch[x]) + '(' + pg(gd_ch, x) + '%)' for x in gd_ch.index]),
                                                             gd_p,
                                                             ';'.join([str(x) + '=' + str(gd_to[x]) + '(' + pg(gd_to, x) + '%)' for x in gd_to.index]),
                                                             ';'.join([str(x) + '=' + str(gd_tn[x]) + '(' + pg(gd_tn, x) + '%)' for x in gd_tn.index]),
                                                             gd_p2,
                                                             'Fisher Exact (2x2)']


# mutation density

md_a = s1_t.groupby('Cohort')['Mutation Density']
md_ch = md_a.get_group('Chapuy et al.').median()
md_sc = md_a.get_group('Schmitz et al.').median()
md_med = s1_t['Mutation Density'].median()

ch_md = s1_t.groupby('Cohort')['Mutation Density'].get_group('Chapuy et al.').dropna()
sc_md = s1_t.groupby('Cohort')['Mutation Density'].get_group('Schmitz et al.').dropna()

md_a = s1_t.groupby('Pair Status (Tumor Only/Normal)')['Mutation Density']
md_tn = md_a.get_group('TN').median()
md_to = md_a.get_group('TO').median()

tn_md = s1_t.groupby('Pair Status (Tumor Only/Normal)')['Mutation Density'].get_group('TN').dropna()
to_md = s1_t.groupby('Pair Status (Tumor Only/Normal)')['Mutation Density'].get_group('TO').dropna()

_, md_p = ss.mannwhitneyu(ch_md, sc_md, method="asymptotic")
_, md_p2 = ss.mannwhitneyu(tn_md, to_md, method="asymptotic")

md_med = '%.3E' % Decimal(np.format_float_positional(md_med, precision=3, fractional=False, min_digits=2))
md_sc = '%.3E' % Decimal(np.format_float_positional(md_sc, precision=3, fractional=False, min_digits=2))
md_ch = '%.3E' % Decimal(np.format_float_positional(md_ch, precision=3, fractional=False, min_digits=2))
md_to = '%.3E' % Decimal(np.format_float_positional(md_to, precision=3, fractional=False, min_digits=2))
md_tn = '%.3E' % Decimal(np.format_float_positional(md_tn, precision=3, fractional=False, min_digits=2))

table_s2b.loc['Mutation Density (Median)'] = [md_med,
                                              md_sc, md_ch, md_p,
                                              md_to, md_tn, md_p2,
                                              'Mann-Whitney U']

# Mutation Density Non Silent
md_a = s1_t.groupby('Cohort')['Mutation Density (Non Silent)']
md_ch = md_a.get_group('Chapuy et al.').median()
md_sc = md_a.get_group('Schmitz et al.').median()
md_med = s1_t['Mutation Density (Non Silent)'].median()

ch_md = s1_t.groupby('Cohort')['Mutation Density (Non Silent)'].get_group('Chapuy et al.').dropna()
sc_md = s1_t.groupby('Cohort')['Mutation Density (Non Silent)'].get_group('Schmitz et al.').dropna()

md_a = s1_t.groupby('Pair Status (Tumor Only/Normal)')['Mutation Density (Non Silent)']
md_tn = md_a.get_group('TN').median()
md_to = md_a.get_group('TO').median()

tn_md = s1_t.groupby('Pair Status (Tumor Only/Normal)')['Mutation Density (Non Silent)'].get_group('TN').dropna()
to_md = s1_t.groupby('Pair Status (Tumor Only/Normal)')['Mutation Density (Non Silent)'].get_group('TO').dropna()

_, md_p = ss.mannwhitneyu(ch_md, sc_md, method="asymptotic")
_, md_p2 = ss.mannwhitneyu(tn_md, to_md, method="asymptotic")

md_med = '%.3E' % Decimal(np.format_float_positional(md_med, precision=3, fractional=False, min_digits=2))
md_sc = '%.3E' % Decimal(np.format_float_positional(md_sc, precision=3, fractional=False, min_digits=2))
md_ch = '%.3E' % Decimal(np.format_float_positional(md_ch, precision=3, fractional=False, min_digits=2))
md_to = '%.3E' % Decimal(np.format_float_positional(md_to, precision=3, fractional=False, min_digits=2))
md_tn = '%.3E' % Decimal(np.format_float_positional(md_tn, precision=3, fractional=False, min_digits=2))


table_s2b.loc['Mutation Density (Non Silent, Median)'] = [md_med,
                                                          md_sc, md_ch, md_p,
                                                          md_to, md_tn, md_p2,
                                                          'Mann-Whitney U']
# driver counts muts
dc_a = s1_t.groupby('Cohort')['Number of Drivers - Mutations']
dc_ch = int(dc_a.get_group('Chapuy et al.').sum())
dc_sc = int(dc_a.get_group('Schmitz et al.').sum())
dc_f = int(s1_t['Number of Drivers - Mutations'].sum())

dc_a = s1_t.groupby('Cohort')['Number of WT Mutations (0)']
wt_ch = int(dc_a.get_group('Chapuy et al.').sum())
wt_sc = int(dc_a.get_group('Schmitz et al.').sum())
wt_f = int(s1_t['Number of WT Mutations (0)'].sum())

sc_ch_t = np.matrix([[dc_ch, wt_ch], [dc_sc, wt_sc]])
dc_p = fisher_exact_2x2(sc_ch_t)[0][0]

dc_a = s1_t.groupby('Pair Status (Tumor Only/Normal)')['Number of Drivers - Mutations']
dc_tn = int(dc_a.get_group('TN').sum())
dc_to = int(dc_a.get_group('TO').sum())

dc_a = s1_t.groupby('Pair Status (Tumor Only/Normal)')['Number of WT Mutations (0)']
wt_tn = int(dc_a.get_group('TN').sum())
wt_to = int(dc_a.get_group('TO').sum())

to_tn_t = np.matrix([[dc_tn, wt_tn], [dc_to, wt_to]])
dc_p2 = fisher_exact_2x2(to_tn_t)[0][0]

table_s2b.loc['Driver Count (Mutations, GSM)'] = ['V=' + str(dc_f) + ';WT=' + str(wt_f),
                                                         'V=' + str(dc_sc) + ';WT=' + str(wt_sc),
                                                         'V=' + str(dc_ch) + ';WT=' + str(wt_ch),
                                                         dc_p,
                                                         'V=' + str(dc_to) + ';WT=' + str(wt_to),
                                                         'V=' + str(dc_tn) + ';WT=' + str(wt_tn),
                                                         dc_p2,
                                                         'Fisher Exact (2x2)']

# driver counts scnas
dc_a = s1_t.groupby('Cohort')['Number of Drivers - SCNAs']
dc_ch = int(dc_a.get_group('Chapuy et al.').sum())
dc_sc = int(dc_a.get_group('Schmitz et al.').sum())
dc_f = int(s1_t['Number of Drivers - SCNAs'].sum())

dc_a = s1_t.groupby('Cohort')['Number of WT SCNAs (0)']
wt_ch = int(dc_a.get_group('Chapuy et al.').sum())
wt_sc = int(dc_a.get_group('Schmitz et al.').sum())
wt_f = int(s1_t['Number of WT SCNAs (0)'].sum())

sc_ch_t = np.matrix([[dc_ch, wt_ch], [dc_sc, wt_sc]])
dc_p = fisher_exact_2x2(sc_ch_t)[0][0]

dc_a = s1_t.groupby('Pair Status (Tumor Only/Normal)')['Number of Drivers - SCNAs']
dc_tn = int(dc_a.get_group('TN').sum())
dc_to = int(dc_a.get_group('TO').sum())

dc_a = s1_t.groupby('Pair Status (Tumor Only/Normal)')['Number of WT SCNAs (0)']
wt_tn = int(dc_a.get_group('TN').sum())
wt_to = int(dc_a.get_group('TO').sum())

to_tn_t = np.matrix([[dc_tn, wt_tn], [dc_to, wt_to]])
dc_p2 = fisher_exact_2x2(to_tn_t)[0][0]

table_s2b.loc['Driver Count (SCNAs, GSM)'] = ['V=' + str(dc_f) + ';WT=' + str(wt_f),
                                              'V=' + str(dc_sc) + ';WT=' + str(wt_sc),
                                              'V=' + str(dc_ch) + ';WT=' + str(wt_ch),
                                              dc_p,
                                              'V=' + str(dc_to) + ';WT=' + str(wt_to),
                                              'V=' + str(dc_tn) + ';WT=' + str(wt_tn),
                                              dc_p2,
                                              'Fisher Exact (2x2)']

# driver counts deletions
dc_a = s1_t.groupby('Cohort')['Number of Drivers - Deletions']
dc_ch = int(dc_a.get_group('Chapuy et al.').sum())
dc_sc = int(dc_a.get_group('Schmitz et al.').sum())
dc_f = int(s1_t['Number of Drivers - Deletions'].sum())

dc_a = s1_t.groupby('Cohort')['Number of WT Deletions (0)']
wt_ch = int(dc_a.get_group('Chapuy et al.').sum())
wt_sc = int(dc_a.get_group('Schmitz et al.').sum())
wt_f = int(s1_t['Number of WT Deletions (0)'].sum())

sc_ch_t = np.matrix([[dc_ch, wt_ch], [dc_sc, wt_sc]])
dc_p = fisher_exact_2x2(sc_ch_t)[0][0]

dc_a = s1_t.groupby('Pair Status (Tumor Only/Normal)')['Number of Drivers - Deletions']
dc_tn = int(dc_a.get_group('TN').sum())
dc_to = int(dc_a.get_group('TO').sum())

dc_a = s1_t.groupby('Pair Status (Tumor Only/Normal)')['Number of WT Deletions (0)']
wt_tn = int(dc_a.get_group('TN').sum())
wt_to = int(dc_a.get_group('TO').sum())

to_tn_t = np.matrix([[dc_tn, wt_tn], [dc_to, wt_to]])
dc_p2 = fisher_exact_2x2(to_tn_t)[0][0]

table_s2b.loc['Driver Count (Deletions, GSM)'] = ['V=' + str(dc_f) + ';WT=' + str(wt_f),
                                                  'V=' + str(dc_sc) + ';WT=' + str(wt_sc),
                                                  'V=' + str(dc_ch) + ';WT=' + str(wt_ch),
                                                  dc_p,
                                                  'V=' + str(dc_to) + ';WT=' + str(wt_to),
                                                  'V=' + str(dc_tn) + ';WT=' + str(wt_tn),
                                                  dc_p2,
                                                  'Fisher Exact (2x2)']

# driver counts amplifications
dc_a = s1_t.groupby('Cohort')['Number of Drivers - Amplifications']
dc_ch = int(dc_a.get_group('Chapuy et al.').sum())
dc_sc = int(dc_a.get_group('Schmitz et al.').sum())
dc_f = int(s1_t['Number of Drivers - Amplifications'].sum())

dc_a = s1_t.groupby('Cohort')['Number of WT Amplifications (0)']
wt_ch = int(dc_a.get_group('Chapuy et al.').sum())
wt_sc = int(dc_a.get_group('Schmitz et al.').sum())
wt_f = int(s1_t['Number of WT Amplifications (0)'].sum())

sc_ch_t = np.matrix([[dc_ch, wt_ch], [dc_sc, wt_sc]])
dc_p = fisher_exact_2x2(sc_ch_t)[0][0]

dc_a = s1_t.groupby('Pair Status (Tumor Only/Normal)')['Number of Drivers - Amplifications']
dc_tn = int(dc_a.get_group('TN').sum())
dc_to = int(dc_a.get_group('TO').sum())

dc_a = s1_t.groupby('Pair Status (Tumor Only/Normal)')['Number of WT Amplifications (0)']
wt_tn = int(dc_a.get_group('TN').sum())
wt_to = int(dc_a.get_group('TO').sum())

to_tn_t = np.matrix([[dc_tn, wt_tn], [dc_to, wt_to]])
dc_p2 = fisher_exact_2x2(to_tn_t)[0][0]

table_s2b.loc['Driver Count (Amplifications, GSM)'] = ['V=' + str(dc_f) + ';WT=' + str(wt_f),
                                                       'V=' + str(dc_sc) + ';WT=' + str(wt_sc),
                                                       'V=' + str(dc_ch) + ';WT=' + str(wt_ch),
                                                       dc_p,
                                                       'V=' + str(dc_to) + ';WT=' + str(wt_to),
                                                       'V=' + str(dc_tn) + ';WT=' + str(wt_tn),
                                                       dc_p2,
                                                       'Fisher Exact (2x2)']

# driver counts scnas
dc_a = s1_t.groupby('Cohort')['Number of Drivers - SVs']
dc_ch = int(dc_a.get_group('Chapuy et al.').sum())
dc_sc = int(dc_a.get_group('Schmitz et al.').sum())
dc_f = int(s1_t['Number of Drivers - SVs'].sum())

dc_a = s1_t.groupby('Cohort')['Number of WT SVs (0)']
wt_ch = int(dc_a.get_group('Chapuy et al.').sum())
wt_sc = int(dc_a.get_group('Schmitz et al.').sum())
wt_f = int(s1_t['Number of WT SVs (0)'].sum())

sc_ch_t = np.matrix([[dc_ch, wt_ch], [dc_sc, wt_sc]])
dc_p = fisher_exact_2x2(sc_ch_t)[0][0]

dc_a = s1_t.groupby('Pair Status (Tumor Only/Normal)')['Number of Drivers - SVs']
dc_tn = int(dc_a.get_group('TN').sum())
dc_to = int(dc_a.get_group('TO').sum())

dc_a = s1_t.groupby('Pair Status (Tumor Only/Normal)')['Number of WT SVs (0)']
wt_tn = int(dc_a.get_group('TN').sum())
wt_to = int(dc_a.get_group('TO').sum())

to_tn_t = np.matrix([[dc_tn, wt_tn], [dc_to, wt_to]])
dc_p2 = fisher_exact_2x2(to_tn_t)[0][0]

table_s2b.loc['Driver Count (SVs, GSM)'] = ['V=' + str(dc_f) + ';WT=' + str(wt_f),
                                            'V=' + str(dc_sc) + ';WT=' + str(wt_sc),
                                            'V=' + str(dc_ch) + ';WT=' + str(wt_ch),
                                            dc_p,
                                            'V=' + str(dc_to) + ';WT=' + str(wt_to),
                                            'V=' + str(dc_tn) + ';WT=' + str(wt_tn),
                                            dc_p2,
                                            'Fisher Exact (2x2)']


table_s2b['p value (Schmitz vs Chapuy)'] = table_s2b['p value (Schmitz vs Chapuy)'].apply(lambda x: x if x == '' else '%.3E' %
                                                                                                  Decimal(np.format_float_positional(x,
                                                                                                                                     precision=3,
                                                                                                                                     fractional=False,
                                                                                                                                     min_digits=2)))

table_s2b['p value (TO vs TN)'] = table_s2b['p value (TO vs TN)'].apply(lambda x: x if x == '' else '%.3E' %
                                                                                         Decimal(np.format_float_positional(x,
                                                                                                                            precision=3,
                                                                                                                            fractional=False,
                                                                                                                            min_digits=2)))

table_s2b.to_csv('../../data_tables/tableS2B_genomic.tsv', sep='\t')
