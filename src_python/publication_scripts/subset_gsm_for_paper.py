import pandas as pd

gsm = pd.read_csv('../../data_tables/gsm/DLBCL_Staudt_Shipp_CL.for_classifier_training.classifier_subset.fix_sv.fix_ploidy.17-Aug-2022.txt',
                  sep='\t', index_col=0)
sig_genes = pd.read_csv('../../data_tables/qval_dfs/fisher_exact_5x2_17-Aug-2022.combined.tsv', sep='\t', index_col=0)
sig_genes = sig_genes.loc[sig_genes['q'] < 0.10].index

gsm_sub = gsm.loc[sig_genes]

gsm_sub = gsm_sub.astype(float).astype(int)
gsm_sub.insert(0, 'classifier_name', gsm_sub.index)

gsm_sub.index = [x[1::] if '.DEL' in x or '.AMP' in x else x for x in gsm_sub.index]

gsm_sub.to_csv('../../data_tables/gsm/GSM_s699_f162_Aug_17_2022.tsv', sep='\t', index_label='gene')