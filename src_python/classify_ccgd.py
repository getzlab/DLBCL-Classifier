import pandas as pd
import numpy as np
import nn
import random
import format_data
import classify_generic as CG

random.seed(123)
np.random.seed(123)

datafile = "../data_tables/gsm/GSM.CCGD.updated.Aug-17-2022.tsv"
df = pd.read_csv(datafile, sep='\t', index_col=0).T

modelpath = '../saved_models/FINALMODEL_NN_evaluation_seeds1_100_folds5_reducedV3.2_removeN5/*'

reduced_df = format_data.construct_reduced_winning_version(df)
# row = {'BCL6_ALT': 0, 'NOTCH2_vec': 1, 'M88O_vec': 0, 'C1_vec4': 0, 'CD70_vec': 0, 'TP53_biallelic':2,
#        'X21Q_AMP': 0, 'Sum_C2_ARM': 4, 'Sum_C2_FOCAL': 12, 'BCL2_combined': 0, 'CREBBP_vec': 0, 'GNA13_vec': 0, 'PTEN': 0, 'SV_MYC': 0,
#        'Hist_comp': 0,	'SGK1_vec': 0, 'ZFP36L1_vec': 4, 'CN_2P16_1_AMP': 0, 'TBL1XR1_vec': 11, 'MYD88_L265P_CD79B': 4, 'Sum_C5_CNA': 6}
#
# hbl1_df = pd.DataFrame.from_dict(row, orient='index').T
# reduced_df = reduced_df.append(hbl1_df)
#
# reduced_df.loc['HBL1_old'] = reduced_df.loc[0]
# reduced_df = reduced_df.drop(0, axis=0)
#
#
# reduced_df.loc['HBL1_notp53'] = reduced_df.loc['HBL1']
# reduced_df.loc['HBL1_notp53', 'TP53_biallelic'] = 0

pred_df = CG.classify_samples_winning_model(reduced_df)

reduced_df.to_csv('../evaluation_cell_lines/ccgd_lines_reducedV3.2_removeN5_gsm.tsv', sep='\t')
pred_df = pred_df.sort_values(by=['PredictedCluster', 'Confidence'])
pred_df.to_csv('../evaluation_cell_lines/CCGD_predictions.tsv', sep='\t')
