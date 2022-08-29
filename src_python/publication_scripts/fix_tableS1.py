import pandas as pd
import numpy as np

table_s1_bjoern = pd.read_csv('../../data_tables/tableS1_bjoern.tsv', sep='\t', index_col=2)
current_s1 = pd.read_csv('../../data_tables/tableS1_classifier.tsv', sep='\t', index_col=0)
staudt_s9 = pd.read_csv('../../data_tables/Tab S9 Characteristics DLBCL-Table 1.tsv', sep='\t', index_col=0)

cols_to_keep = ['OriginalCluster [chapuy et al Nature Med 2018]', 'Above85', 'Above80',
                'order', 'Witin.700', 'MatchID', 'Gender', 'Age-at first diagnosis',
                'IPI_AGE', 'IPI_ECOG', 'IPI_STAGE', 'IPI_LDH',
                'IPI_EXBM', 'IPI', 'Biopsy Type', 'R-CHOp-like\nChemo', 'PFS-years', 'PFS_STAT',
                'OS.time', 'OS.status (1=dead)', 'LymphGenClass [Schmitz]\ncall', 'LymphGenClass [Wright]\ncall',
                'Dbl.Hit\nCall', 'ob24z (CNS involvment)', 'hodz (Testicular invovlement)']


DLBCL11454 = 'F'
DLBCL11486 = 'M'
DLBCL11581 = 'M'

current_s1 = current_s1.loc[table_s1_bjoern.index, ~current_s1.columns.isin(cols_to_keep)]

combined_s1 = pd.concat([current_s1, table_s1_bjoern[cols_to_keep]], axis=1)
combined_s1 = combined_s1.replace('^na$', np.nan, regex=True)

combined_s1.loc['DLBCL11454', 'Gender'] = DLBCL11454
combined_s1.loc['DLBCL11486', 'Gender'] = DLBCL11486
combined_s1.loc['DLBCL11581', 'Gender'] = DLBCL11581

combined_s1 = combined_s1.rename({'Cluster': 'Cluster [NMF new]'}, axis=1)
combined_s1.index.name = 'ID'

combined_s1.isna().sum(axis=0).sort_values().to_csv('../../data_tables/na_counts_tables1.tsv', sep='\t')

# Index(['Cohort', 'COO', 'Sample Preparation', 'Mean Target Coverage',
#        'Median Target Coverage', 'Ploidy (ABSOLUTE)', 'Purity (ABSOLUTE)',
#        'Pair Status (Tumor Only/Normal)', 'Number of Mutations',
#        'Fraction Genome Deleted', 'Fraction Genome Amplified',
#        'Number of Drivers - Mutations',
#        'Number of Drivers - Non-Silent Mutations', 'Number of Drivers - SCNAs',
#        'Test/Train Set Membership', 'Cluster [NMF new]', 'Target C1',
#        'Target C2', 'Target C3', 'Target C4', 'Target C5', 'Predicted C1',
#        'Predicted C2', 'Predicted C3', 'Predicted C4', 'Predicted C5',
#        'PredictedCluster', 'Confidence', 'Above90', 'Top70 Perc. Confident',
#        'Staudt Only Cluster (k=5)', 'Staudt Only Cluster (k=4)',
#        'Shipp Only Cluster', 'OriginalCluster [chapuy et al Nature Med 2018]',
#        'Above85', 'Above80',
#        'Age-at first diagnosis', 'IPI_AGE', 'IPI_ECOG', 'IPI_STAGE', 'IPI_LDH',
#        'IPI_EXBM', 'IPI', 'Biopsy Type', 'R-CHOp-like\nChemo', 'PFS-years',
#        'PFS_STAT', 'OS.time', 'OS.status (1=dead)',
#        'LymphGenClass [Schmitz]\ncall', 'LymphGenClass [Wright]\ncall',
#        'Dbl.Hit\nCall', 'ob24z (CNS involvment)',
#        'hodz (Testicular invovlement)'],
#       dtype='object')

# add mut density

mut_dens = pd.read_csv('../../data_tables/mutsig2cv_gistic_qvalues/DLBCL_551_training.patient_counts_and_rates.txt', sep='\t', index_col=0)
mut_dens = mut_dens.drop('DLBCL11493')
combined_s1['Mutation Density'] = np.nan
combined_s1['Mutation Density (Non Silent)'] = np.nan

combined_s1.loc[mut_dens.index, 'Mutation Density'] = mut_dens['rate_tot']
combined_s1.loc[mut_dens.index, 'Mutation Density (Non Silent)'] = mut_dens['rate_non']


col_order = ['order', 'MatchID', 'Gender', 'Cohort', 'COO', 'Sample Preparation', 'Mean Target Coverage',
             'Median Target Coverage', 'Ploidy (ABSOLUTE)', 'Purity (ABSOLUTE)', 'Pair Status (Tumor Only/Normal)',
             'Number of Mutations', 'Fraction Genome Deleted', 'Fraction Genome Amplified',
             'Number of Drivers - Mutations', 'Number of Drivers - Non-Silent Mutations',
             'Number of Drivers - SCNAs', 'Number of Drivers - Amplifications', 'Number of Drivers - Deletions',
             'Number of Drivers - SVs',
             'Number of WT Mutations (0)', 'Number of WT SCNAs (0)', 'Number of WT Amplifications (0)', 'Number of WT Deletions (0)',
             'Number of WT SVs (0)',
             'Mutation Density', 'Mutation Density (Non Silent)',
             'Age-at first diagnosis', 'IPI_AGE', 'IPI_ECOG', 'IPI_STAGE', 'IPI_LDH',
             'IPI_EXBM', 'IPI', 'Biopsy Type', 'R-CHOp-like\nChemo', 'PFS-years',
             'PFS_STAT', 'OS.time', 'OS.status (1=dead)','LymphGenClass [Schmitz]\ncall', 'LymphGenClass [Wright]\ncall',
             'Dbl.Hit\nCall', 'ob24z (CNS involvment)',
             'hodz (Testicular invovlement)',
             'Test/Train Set Membership', 'Cluster [NMF new]', 'Target C1',
             'Target C2', 'Target C3', 'Target C4', 'Target C5', 'Predicted C1',
             'Predicted C2', 'Predicted C3', 'Predicted C4', 'Predicted C5',
             'PredictedCluster', 'Confidence', 'Above90', 'Above85', 'Above80',
             'Top70 Perc. Confident', 'Staudt Only Cluster',
             'Shipp Only Cluster', 'OriginalCluster [chapuy et al Nature Med 2018]'
             ]

staudt_s9 = staudt_s9.loc[staudt_s9.index.isin(combined_s1.index)]
# fix IPI

# apply fixes to weird values (fractional ECOG, -1 LDH ratio, etc).
ecog_vals = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0}

staudt_s9.loc[staudt_s9['LDH Ratio'] < 0, 'LDH Ratio'] = np.nan
staudt_s9.loc[~staudt_s9['ECOG Performance Status'].isin(ecog_vals), 'ECOG Performance Status'] = np.nan
staudt_s9.loc[staudt_s9['Number of Extranodal Sites'] % 1 != 0, 'Number of Extranodal Sites'] = np.nan

age_bool = staudt_s9['Age'].apply(lambda x: pd.NA if np.isnan(x) else int(x >= 60))
aa_bool = staudt_s9['Ann Arbor Stage'].apply(lambda x: pd.NA if np.isnan(x) else int(x > 2))
ldh_bool = staudt_s9['LDH Ratio'].apply(lambda x: pd.NA if np.isnan(x) else int(x > 2))
ecog_bool = staudt_s9['ECOG Performance Status'].apply(lambda x: pd.NA if np.isnan(x) else int(x > 2))
en_bool = staudt_s9['Number of Extranodal Sites'].apply(lambda x: pd.NA if np.isnan(x) else int(x > 2))


combined_s1.loc[age_bool.index, 'IPI'] = age_bool + aa_bool + ldh_bool + ecog_bool + en_bool
combined_s1.loc[age_bool.index, 'IPI_AGE'] = age_bool
combined_s1.loc[ecog_bool.index, 'IPI_ECOG'] = ecog_bool
combined_s1.loc[aa_bool.index, 'IPI_STAGE'] = aa_bool
combined_s1.loc[ldh_bool.index, 'IPI_LDH'] = ldh_bool
combined_s1.loc[en_bool.index, 'IPI_EXBM'] = en_bool


combined_s1 = combined_s1[col_order]
combined_s1.to_csv('../../data_tables/tableS1_classifier_merged.tsv', sep='\t')
