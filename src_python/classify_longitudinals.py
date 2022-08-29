import pandas as pd
import numpy as np
import nn
import random
import format_data
import classify_generic as CG

random.seed(123)
np.random.seed(123)

############
# Classify #
############

datafile = '../data_tables/gsm/Quality_Longitudinals.01-Sep-2021.txt'
qval_file = '../data_tables/qval_dfs/fisher_exact_5x2_17-Aug-2022.combined.tsv'
df = pd.read_csv(datafile, sep='\t', index_col=0).T
df = df.drop('cohort')
df.columns = df.columns.str.replace('_', '.')
df.columns = df.columns.str.replace('-', '.')
df.columns = ['X' + x if '.AMP' in x or '.DEL' in x else x for x in df.columns]
df['MYD88'] = df[['MYD88.L265P', 'MYD88.OTHER']].max(axis=1)

reduced_df = format_data.construct_reduced_winning_version(df)
pred_df = CG.classify_samples_winning_model(reduced_df)

##########################################
# Add longitudinal metadata to the table #
##########################################


datafile = '../data_tables/purities_ploidies/Quality_Longitudinals.pairs.txt'
lng_meta = pd.read_csv(datafile, sep='\t', index_col=0)
lng_meta.index = lng_meta.index.str.replace('-', '_')

pred_df['purity'] = lng_meta['purity_absolute_reviewed']
pred_df['ploidy'] = df['PLOIDY']

drivers = pd.read_csv(qval_file, sep='\t', index_col=0)

df_drivers = df[drivers.loc[drivers['q'] <= 0.10].index]
df_classifier = df[drivers.loc[drivers['q'] <= 0.05].index]

sum_all = (df_drivers != 0).sum(axis=1)
sum_classifier = (df_classifier != 0).sum(axis=1)

pred_df['sum_drivers_all'] = sum_all
pred_df['sum_drivers_classifier'] = sum_classifier

pred_df = pd.concat([pred_df, reduced_df], axis=1)

pred_df = pred_df.loc[pred_df['purity'] >= 0.20]

###############
# Write files #
###############

reduced_df.to_csv('../evaluation_longitudinals/longitudinals_reduced_gsm.tsv', sep='\t')
pred_df = pred_df.sort_values(by=['PredictedCluster', 'Confidence'])
pred_df.to_csv('../evaluation_longitudinals/longitudinals_predictions.tsv', sep='\t')
pred_df_sortedname = pred_df.sort_index()
pred_df_sortedname.to_csv('../evaluation_longitudinals/longitudinals_predictions_namesorted.tsv', sep='\t')
