import pandas as pd

all_samples = pd.read_csv('../../data_tables/sample_sets/ShippStaudtSets.purity0.2.txt',
                          sep='\t', index_col=0)
shipp_samples = all_samples.loc[all_samples['cohort'] == 'Shipp'].index
staudt_samples = all_samples.loc[all_samples['cohort'] != 'Shipp'].index

full_gsm = pd.read_csv('../../data_tables/gsm/DLBCL.699.fullGSM.Sep_23_2022.tsv', sep='\t', index_col=0)
nm_gsm = pd.read_csv('../../data_tables/gsm/DLBCL_significant_event_matrix_NatMed.txt',
                     sep='\t', index_col=0)
nm_gsm.index = nm_gsm.index.str.replace(':', '.').str.replace('-', '.').str.upper()
nm_gsm.index = ['X' + x if '.AMP' in x or '.DEL' in x else x for x in nm_gsm.index]
nm_gsm.index = [x if x != 'X18Q21.33(BCL2).AMP' else 'X18Q21.32.AMP' for x in nm_gsm.index]

clus_gsm = full_gsm.loc[full_gsm.index.isin(nm_gsm.index)]
clus_gsm_shipp = clus_gsm.loc[:, clus_gsm.columns[clus_gsm.columns.isin(shipp_samples)]]
clus_gsm_staudt = clus_gsm.loc[:, clus_gsm.columns[clus_gsm.columns.isin(staudt_samples)]]

clus_gsm.index.name = 'gene'
clus_gsm_shipp.index.name = 'gene'
clus_gsm_staudt.index.name = 'gene'

clus_gsm.to_csv('../../data_tables/gsm/DLBCL.699.fullGSM.Sep_23_2022.CLUSTERING.tsv',
                sep='\t')
clus_gsm_shipp.to_csv('../../data_tables/gsm/DLBCL.' + str(clus_gsm_shipp.shape[1]) + '.fullGSM.Sep_23_2022.CLUSTERING.SHIPP.tsv',
                      sep='\t')
clus_gsm_staudt.to_csv('../../data_tables/gsm/DLBCL.' + str(clus_gsm_staudt.shape[1]) + '.fullGSM.Sep_23_2022.CLUSTERING.STAUDT.tsv',
                       sep='\t')
