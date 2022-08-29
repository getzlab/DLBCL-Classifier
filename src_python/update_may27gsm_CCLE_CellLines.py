import pandas as pd
import re


GSM = pd.read_csv('../data_tables/gsm/GSM.depmap.lymphoma.20Jun2020.tsv', sep='\t', index_col=0)
GSM.drop('Driver_Discovery', axis=1, inplace=True)

dlbcl_mat = pd.read_csv("../data_tables/gsm/DLBCL_Staudt_Shipp_CL.for_classifier_training.classifier_subset.fix_sv.fix_ploidy.17-Aug-2022.txt",
                        sep='\t', index_col=0)

qval_file = '../data_tables/qval_dfs/fisher_exact_5x2_17-Aug-2022.combined.tsv'
qval_df = pd.read_csv(qval_file, sep='\t', index_col=0)
qval_df = qval_df.loc[qval_df['q'] <= 0.10]

depmapMAF = pd.read_csv('../data_tables/maf_files/CCLE_mutations_trimmed.csv', sep='\t', low_memory=False)
samplemapping = pd.read_csv('../data_tables/maf_files/sample_info.csv', sep=',')
samplemapping = pd.Series(samplemapping['CCLE_Name'].values, index=samplemapping['DepMap_ID']).to_dict()

sv_indels = [x for x in dlbcl_mat.index if '.AMP' in x or '.DEL' in x or 'SV.' in x]
GSM = GSM.loc[GSM.index.isin(sv_indels)]

genes_to_add = [x for x in qval_df.index if x not in GSM.index]
genes_to_add = [x.replace('.', '-') for x in genes_to_add]

for gene in genes_to_add:
   GSM.loc[gene] = 0

noncoding_types = {'Intron', "5'UTR", "3'UTR", "5'Flank", "3'Flank", 'IGR', 'COULD_NOT_DETERMINE', 'RNA'}

for idx, row in depmapMAF.iterrows():
    gene = row['Hugo_Symbol']
    if gene not in genes_to_add:
        continue

    if row['Variant_Classification'] in noncoding_types:
        continue

    vc = row['Variant_Classification']

    mappedSample = row['Tumor_Sample_Barcode'].upper()
    if mappedSample in samplemapping:
        mappedSample = samplemapping[mappedSample]

    if mappedSample in GSM.columns:
        if gene == 'MYD88':
            if row['Start_position'] == 38182641 and row['Tumor_Seq_Allele1'] == 'C':
                GSM.loc['MYD88-L265P', mappedSample] = 2
            else:
                if vc == 'Silent':
                    if GSM.loc['MYD88-OTHER', mappedSample] == 0:
                        GSM.loc['MYD88-OTHER', mappedSample] = 1
                else:
                    GSM.loc['MYD88-OTHER', mappedSample] = 2
        else:
            if vc == 'Silent':
                val = 1
            else:
                val = 2
            GSM.at[gene, mappedSample] = max(val, GSM.loc[gene, mappedSample])


depmapSVs = pd.read_csv('../data_tables/gsm/SVs DepMap cell lines.publicdata 4-9-2020.BC.txt', sep='\t', index_col=0)
depmapSVs.index.name = 'Sample'

GSM = GSM.copy(deep=True)
for idx, row in depmapSVs.iterrows():
    GSM.loc['SV.BCL2', idx] = row[1]
    GSM.loc['SV.BCL6', idx] = row[2]
    GSM.loc['SV.MYC', idx] = row[3]

GSM.index = [x.replace('-', '.') for x in GSM.index]
GSM.loc['MYD88'] = GSM.loc[['MYD88.L265P', 'MYD88.OTHER']].max(axis=0)

print(GSM.sum(axis=0).sort_values())
GSM.to_csv('../data_tables/gsm/GSM.Depmap.updated.Aug-17-2022.tsv', sep='\t')

