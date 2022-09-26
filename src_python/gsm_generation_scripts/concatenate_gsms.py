import pandas as pd

OUTPUT_FN = '../../data_tables/gsm/DLBCL.699.fullGSM.Sep_23_2022.tsv'

scna_df = pd.read_csv('../../data_tables/gsm/DLBCL.699.scnaGSM.Sep_23_2022.tsv', sep='\t', index_col=0)
mut_df = pd.read_csv('../../data_tables/gsm/DLBCL.699.mutationsGSM.Sep_23_2022.tsv', sep='\t', index_col=0)
sv_df = pd.read_csv('../../data_tables/gsm/DLBCL.699.svGSM.Sep_23_2022.tsv', sep='\t', index_col=0)

samples = set(scna_df.columns).intersection(set(sv_df.columns)).intersection(set(mut_df.columns))

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

full_gsm = full_gsm[samples]
full_gsm.to_csv(OUTPUT_FN, sep='\t')
