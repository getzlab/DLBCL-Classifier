import pandas as pd


qval_file = '../data_tables/qval_dfs/fisher_exact_5x2.Sep_23_2022.combined.tsv'
longitudinal_gsm_file = "../data_tables/gsm/Quality_Longitudinals.01-Sep-2021.txt"
longitudinal_maf_file = "../data_tables/maf_files/Quality_Longitudinals.maf.fixed_ids.maf"

lng_gsm = pd.read_csv(longitudinal_gsm_file, sep='\t', index_col=0)
lng_maf = pd.read_csv(longitudinal_maf_file, sep='\t')
qval_df = pd.read_csv(qval_file, sep='\t', index_col=0)
sig_genes = qval_df.loc[qval_df['q'] <= 0.10].index

lng_gsm.columns = lng_gsm.columns.str.upper()

lng_gsm = lng_gsm.drop(['COHORT'], axis=1)
lng_gsm = lng_gsm.drop(['EBV'], axis=0)
lng_gsm.index = ['X' + x if '.AMP' in x or '.DEL' in x else x for x in lng_gsm.index]
lng_gsm.index = lng_gsm.index.str.replace('_', '.')
lng_gsm.index = lng_gsm.index.str.replace('-', '.')
lng_gsm = lng_gsm.astype(float).astype(int)

lng_gsm.loc['MYD88'] = lng_gsm.loc[['MYD88.L265P', 'MYD88.OTHER']].max(axis=0)

genes_to_add = [x for x in sig_genes if x not in lng_gsm.index]

for g in genes_to_add:
    lng_gsm.loc[g] = 0

event_subset_maf = lng_maf.loc[lng_maf['Hugo_Symbol'].isin(genes_to_add)]

noncoding_types = {'Intron', "5'UTR", "3'UTR", "5'Flank", "3'Flank", 'IGR', 'COULD_NOT_DETERMINE', 'RNA'}

for i, row in event_subset_maf.iterrows():
    sample = row['Tumor_Sample_Barcode']
    gene = row['Hugo_Symbol']
    vc = row['Variant_Classification']

    if vc in noncoding_types:
        continue

    if vc == 'Silent' and lng_gsm.loc[gene, sample] == 0:
        lng_gsm.loc[gene, sample] = 1
    else:
        lng_gsm.loc[gene, sample] = 2

lng_gsm.to_csv('../data_tables/gsm/Quality_Longitudinals.Sep_23_2022.txt', sep='\t')
