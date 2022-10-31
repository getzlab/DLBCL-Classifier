import pandas as pd
import numpy as np
import nn
import random
import format_data
import classify_generic as CG

random.seed(123)
np.random.seed(123)

datafile = "../data_tables/gsm/GSM.CCGD.updated.Sep_23_2022.tsv"
df = pd.read_csv(datafile, sep='\t', index_col=0).T

reduced_df = format_data.construct_reduced_winning_version(df)

pred_df = CG.classify_samples_winning_model(reduced_df)

reduced_df.to_csv('../evaluation_cell_lines/CCGD_reducedV3.4_removeN5_gsm.tsv', sep='\t')
pred_df = pred_df.sort_values(by=['PredictedCluster', 'Confidence'])
pred_df.to_csv('../evaluation_cell_lines/CCGD_predictions.tsv', sep='\t')
