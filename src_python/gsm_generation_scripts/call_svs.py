import pandas as pd

OUTPUT_FN = '../../data_tables/gsm/DLBCL.699.svGSM.Sep_23_2022.tsv'

sv_file = '../../data_tables/additional_gsm_inputs/DLBCL_Shipp_Staudt.SVs.14-Dec-2021.txt'
sample_set = list(pd.read_csv('../../data_tables/train_test_sets/TrainingSet_550Subset_May2021.txt', sep='\t', header=None, index_col=0).index)
sample_set = set(sample_set + list(pd.read_csv('../../data_tables/train_test_sets/TestingSet_149Subset_May2021.txt', sep='\t', header=None, index_col=0).index))

SV = pd.read_csv(sv_file, sep='\t', index_col=0)
SV['genes'] = SV['gene1'].fillna('') + '>' + SV['gene2'].fillna('')

SV['evt'] = '---'

SV.loc[SV['genes'].str.contains('BCL2'), 'evt'] = 'SV.BCL2'
SV.loc[SV['genes'].str.contains('BCL6'), 'evt'] = 'SV.BCL6'
SV.loc[SV['genes'].str.contains('MYC'), 'evt'] = 'SV.MYC'

SV = SV.loc[SV['evt'] != '---']
SV['CCF'] = SV['CCF'].fillna(1.0)

sv_df = pd.DataFrame(0, columns=sorted(list(sample_set)), index=['SV.BCL2', 'SV.BCL6', 'SV.MYC', 'SV.BCL2.CCF', 'SV.BCL6.CCF', 'SV.MYC.CCF'])
alt_counts_df = pd.DataFrame(0, columns=sorted(list(sample_set)), index=['SV.BCL2', 'SV.BCL6', 'SV.MYC'])

sv_bcl2_samples = set(SV.loc[SV['genes'].str.contains('BCL2'), 'individual']).intersection(sample_set)
sv_myc_samples = set(SV.loc[SV['genes'].str.contains('MYC'), 'individual']).intersection(sample_set)
sv_bcl6_samples = set(SV.loc[SV['genes'].str.contains('BCL6'), 'individual']).intersection(sample_set)

sv_myc_ccfs = SV.loc[SV['individual'].isin(sv_myc_samples) & SV['genes'].str.contains('MYC')].groupby('individual')['CCF'].max()
sv_bcl2_ccfs = SV.loc[SV['individual'].isin(sv_bcl2_samples) & SV['genes'].str.contains('BCL2')].groupby('individual')['CCF'].max()
sv_bcl6_ccfs = SV.loc[SV['individual'].isin(sv_bcl6_samples) & SV['genes'].str.contains('BCL6')].groupby('individual')['CCF'].max()

sv_df.loc['SV.BCL2', sv_bcl2_samples] = 3
sv_df.loc['SV.BCL6', sv_bcl6_samples] = 3
sv_df.loc['SV.MYC', sv_myc_samples] = 3
sv_df.loc['SV.BCL2.CCF', sv_bcl2_ccfs.index] = sv_bcl2_ccfs
sv_df.loc['SV.BCL6.CCF', sv_bcl6_ccfs.index] = sv_bcl6_ccfs
sv_df.loc['SV.MYC.CCF', sv_myc_ccfs.index] = sv_myc_ccfs

for s in sv_bcl2_samples:
    vcf_talt = SV.loc[(SV['individual'] == s) &
                      (SV['genes'].str.contains('BCL2')) &
                      (SV['CCF'] == sv_df.loc['SV.BCL2.CCF', s]), 'VCF_TALT']
    if not vcf_talt.shape[0]:
        print(s)
    alt_counts_df.loc['SV.BCL2', s] = vcf_talt.values[0]

for s in sv_bcl6_samples:
    vcf_talt = SV.loc[(SV['individual'] == s) &
                      (SV['genes'].str.contains('BCL6')) &
                      (SV['CCF'] == sv_df.loc['SV.BCL6.CCF', s]), 'VCF_TALT']
    if not vcf_talt.shape[0]:
        print(s)
    alt_counts_df.loc['SV.BCL6', s] = vcf_talt.values[0]

for s in sv_myc_samples:
    vcf_talt = SV.loc[(SV['individual'] == s) &
                      (SV['genes'].str.contains('MYC')) &
                      (SV['CCF'] == sv_df.loc['SV.MYC.CCF', s]), 'VCF_TALT']
    if not vcf_talt.shape[0]:
        print(s)
    alt_counts_df.loc['SV.MYC', s] = vcf_talt.values[0]

# TCGA SV fix............
tcga_svs = pd.read_csv('../../data_tables/gsm/tcga_svs.tsv', sep='\t', index_col=0)

sv_df[tcga_svs.columns] = tcga_svs
sv_df.to_csv(OUTPUT_FN, sep='\t')
alt_counts_df.to_csv('../../reduce_purity_experiment/sv_alt_counts.tsv', sep='\t')
