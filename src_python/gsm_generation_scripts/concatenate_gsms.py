import pandas as pd

OUTPUT_FN = '../../data_tables/gsm/DLBCL.699.fullGSM.Sep_23_2022.tsv'

scna_df = pd.read_csv('../../data_tables/gsm/DLBCL.699.scnaGSM.Sep_23_2022.tsv', sep='\t', index_col=0)
mut_df = pd.read_csv('../../data_tables/gsm/DLBCL.699.mutationsGSM.Sep_23_2022.tsv', sep='\t', index_col=0)
sv_df = pd.read_csv('../../data_tables/gsm/DLBCL.699.svGSM.Sep_23_2022.tsv', sep='\t', index_col=0)

full_gsm = pd.concat([mut_df, scna_df, sv_df])

# Fix a naming issue
full_gsm.loc['19Q13.32.DEL'] = full_gsm.loc['19Q13.32_1.DEL']
full_gsm.loc['19Q13.32.DEL.CCF'] = full_gsm.loc['19Q13.32_1.DEL.CCF']

full_gsm = full_gsm.drop('19Q13.32_1.DEL', axis=0)
full_gsm = full_gsm.drop('19Q13.32_1.DEL.CCF', axis=0)

# Make names compatible with R... unfortunately..
full_gsm.index = full_gsm.index.str.replace('-', '.')
full_gsm.index = ['X' + x if '.AMP' in x or '.DEL' in x else x for x in full_gsm.index]
full_gsm.loc['COO'] = full_gsm.loc['COO'].str.replace('Unclassified', 'UNC').str.replace('Unclass', 'UNC')

# Set sample & feature order to ensure reproducibility.
# Had to add this after analysis because I was using set() which is not deterministically ordered
sample_order = pd.read_csv('../../data_tables/gsm/sample_order.tsv', index_col=0, header=None).index
feature_order = pd.read_csv('../../data_tables/gsm/feature_order.tsv', index_col=0, header=None).index

full_gsm = full_gsm.loc[feature_order, sample_order]
full_gsm.to_csv(OUTPUT_FN, sep='\t')
