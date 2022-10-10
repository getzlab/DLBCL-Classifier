import argparse
import pandas as pd
import numpy as np

# --maf ../../data_tables/maf_files/DLBCL_combined_699.hg38B.noPDE4DIP.noHISTartifacts.maf --mutsig_sig_genes ../../data_tables/mutsig2cv_gistic_qvalues/DLBCL_550_training_noPDE4DIP_noHISTartifacts.sig_genes.txt --mutsig_q_thresh 0.10 --output_fn ../../data_tables/gsm/DLBCL.699.mutationsGSM.Sep_23_2022.tsv --additional_sig_genes ../../data_tables/additional_gsm_inputs/NatMed_104_sig_genes.tableS3.tsv --include_myd88_L265P --include_ccf CCF_hat --blacklist ../../data_tables/additional_gsm_inputs/dlbclass_blacklist.tsv --ploidy ../../data_tables/purities_ploidies/PloidyDataFrame.txt --purity ../../data_tables/purities_ploidies/ALLPurities_fixednames.tsv --coo ../../data_tables/phenotypes/COO_and_genetic_lables.txt

# this list is hard-defined for now to allow subsetting to non-silent events.
NOT_NONSILENT = ["Silent", "3'UTR", "5'Flank", "5'UTR",
                 "IGR", "Intron", "RNA", "lincRNA",
                 "COULD_NOT_DETERMINE"]

parser = argparse.ArgumentParser()
parser.add_argument('--maf',
                    help='The maf file (TSV format) to call your samples\' mutation events.',
                    required=True, type=str)
parser.add_argument('--mutsig_sig_genes',
                    help='Mutsig sig genes list. GSM events will be called from the significant genes in this list.',
                    required=True, type=str)
parser.add_argument('--mutsig_q_thresh',
                    help='q value to subset significant genes by.',
                    required=True, type=float)
parser.add_argument('--output_fn',
                    help='Output filename.',
                    required=True, type=str)
parser.add_argument('--alt_count_fn',
                    help='Alt counts output filename.',
                    required=False, type=str, default=False)
parser.add_argument('--sample_set',
                    help='A file with one column with the samples to call events for. No header.' +
                         '\nIf none provided, all samples from the column \'Tumor_Sample_Barcode\' will be used',
                    required=False, type=str, default=False)
parser.add_argument('--additional_sig_genes',
                    help='A file where the only column is a list of genes to also include in the sig gene list. No header.',
                    required=False, type=str, default=False)
parser.add_argument('--include_ccf',
                    help='Indicates which column from the maf to include CCF rows for each gene. Does not include CCF if not given.',
                    required=False, type=str, default=False)
parser.add_argument('--include_myd88_L265P',
                    help='Indicates whether or not to call MYD88 as L265P and non-L265P.',
                    required=False, action='store_true', default=False)
parser.add_argument('--blacklist',
                    help='A file where the only column is a list of genes to blacklist (not call). No header.',
                    required=False, type=str, default=False)
parser.add_argument('--ploidy',
                    help='A two column TSV file, with the first column indicating sample name and the second indicating ploidy. No header.',
                    required=False, type=str, default=False)
parser.add_argument('--purity',
                    help='A two column TSV file, with the first column indicating sample name and the second indicating purity. No header.',
                    required=False, type=str, default=False)
parser.add_argument('--coo',
                    help='A TSV file, with two columns named "Sample" and "COO"',
                    required=False, type=str, default=False)

args = parser.parse_args()

cols = ['Hugo_Symbol', 'Tumor_Sample_Barcode', 'Variant_Classification', 'Protein_Change']
if args.include_ccf:
    ccf_col = args.include_ccf
    cols = cols + [ccf_col]

if args.alt_count_fn:
    cols = cols + ['t_alt_count']

maf = pd.read_csv(args.maf, sep='\t', dtype=str, low_memory=False, usecols=cols)

if args.include_ccf:
    maf[ccf_col] = maf[ccf_col].astype(float)

sig_genes = pd.read_csv(args.mutsig_sig_genes, sep='\t', index_col=0)
mutsig_q_thresh = args.mutsig_q_thresh

# Union the sig genes with any additional genes
genes_to_call = set(sig_genes.loc[sig_genes['q'] <= mutsig_q_thresh].index)
additional_genes = set([])
if args.additional_sig_genes:
    additional_genes = set(pd.read_csv(args.additional_sig_genes, header=None).iloc[:, 0].values)

genes_to_call = genes_to_call.union(additional_genes)

# determine which samples to call events for
samples = set(maf['Tumor_Sample_Barcode'])
if args.sample_set:
    samples = sorted(list(set(pd.read_csv(args.sample_set, header=None).iloc[:, 0].values)))

# determine blacklist
blacklist = set([])
if args.blacklist:
    blacklist = set(pd.read_csv(args.blacklist, header=None).iloc[:, 0].values)

genes_to_call = genes_to_call - blacklist

# make a samples x genes_to_call GSM
GSM = pd.DataFrame(0, index=genes_to_call, columns=samples)

if args.include_myd88_L265P:
    GSM.loc['MYD88.L265P'] = 0
    GSM.loc['MYD88.OTHER'] = 0

if args.alt_count_fn:
    alt_count_df = GSM.copy(deep=True)

if args.include_ccf:
    ccf_GSM = pd.DataFrame(-1, index=GSM.index, columns=samples)
    ccf_GSM.index = ccf_GSM.index + '.CCF'
    GSM = pd.concat([GSM, ccf_GSM])

# Fill in non-silent events first.

non_sil_maf = maf.loc[(maf['Hugo_Symbol'].isin(genes_to_call)) &
                      ~(maf['Variant_Classification'].isin(NOT_NONSILENT))]

for _, event in non_sil_maf.iterrows():
    gene = event['Hugo_Symbol']
    s = event['Tumor_Sample_Barcode']

    GSM.loc[gene, s] = 2

    if args.include_myd88_L265P and gene == 'MYD88':
        if event['Protein_Change'] == 'p.L265P':
            GSM.loc['MYD88.L265P', s] = 2
            if args.include_ccf:
                if np.isnan(GSM.loc['MYD88.L265P' + '.CCF', s]) or GSM.loc['MYD88.L265P' + '.CCF', s] == -1:
                    GSM.loc['MYD88.L265P' + '.CCF', s] = event[ccf_col]
                    if args.alt_count_fn:
                        alt_count_df.loc['MYD88.L265P', s] = event['t_alt_count']
                elif not np.isnan(event[ccf_col]):
                    GSM.loc['MYD88.L265P' + '.CCF', s] = max(GSM.loc['MYD88.L265P' + '.CCF', s], event[ccf_col])
                    if event[ccf_col] > GSM.loc['MYD88.L265P' + '.CCF', s]:
                        alt_count_df.loc['MYD88.L265P', s] = event['t_alt_count']
        else:
            GSM.loc['MYD88.OTHER', s] = 2
            if args.include_ccf:
                if np.isnan(GSM.loc['MYD88.OTHER' + '.CCF', s]) or GSM.loc['MYD88.OTHER' + '.CCF', s] == -1:
                    GSM.loc['MYD88.OTHER' + '.CCF', s] = event[ccf_col]
                    if args.alt_count_fn:
                        alt_count_df.loc['MYD88.OTHER', s] = event['t_alt_count']
                elif not np.isnan(event[ccf_col]):
                    GSM.loc['MYD88.OTHER' + '.CCF', s] = max(GSM.loc['MYD88.OTHER' + '.CCF', s], event[ccf_col])
                    if event[ccf_col] > GSM.loc['MYD88.OTHER' + '.CCF', s]:
                        alt_count_df.loc['MYD88.OTHER', s] = event['t_alt_count']

    if args.include_ccf:
        if np.isnan(GSM.loc[gene + '.CCF', s]) or GSM.loc[gene + '.CCF', s] == -1:
            GSM.loc[gene + '.CCF', s] = event[ccf_col]
            if args.alt_count_fn:
                alt_count_df.loc[gene, s] = event['t_alt_count']
        elif not np.isnan(event[ccf_col]):
            GSM.loc[gene + '.CCF', s] = max(GSM.loc[gene + '.CCF', s], event[ccf_col])
            if event[ccf_col] > GSM.loc[gene + '.CCF', s]:
                alt_count_df.loc[gene, s] = event['t_alt_count']


            # Next call silents. Don't clobber non silent events if there.

sil_maf = maf.loc[(maf['Hugo_Symbol'].isin(genes_to_call)) &
                  (maf['Variant_Classification'] == 'Silent')]

for _, event in sil_maf.iterrows():
    gene = event['Hugo_Symbol']
    s = event['Tumor_Sample_Barcode']

    if gene == 'MYD88':
        GSM.loc['MYD88.OTHER', s] = 1
        if args.include_ccf:
            if np.isnan(GSM.loc['MYD88.OTHER' + '.CCF', s]) or GSM.loc['MYD88.OTHER' + '.CCF', s] == -1:
                GSM.loc['MYD88.OTHER' + '.CCF', s] = event[ccf_col]
                if args.alt_count_fn:
                    alt_count_df.loc['MYD88.OTHER', s] = event['t_alt_count']
            elif not np.isnan(event[ccf_col]):
                GSM.loc['MYD88.OTHER' + '.CCF', s] = max(GSM.loc['MYD88.OTHER' + '.CCF', s], event[ccf_col])
                if event[ccf_col] > GSM.loc['MYD88.OTHER' + '.CCF', s]:
                    alt_count_df.loc['MYD88.OTHER', s] = event['t_alt_count']

    if GSM.loc[gene, s] == 2:
        continue

    GSM.loc[gene, s] = 1

    if args.include_ccf:
        if np.isnan(GSM.loc[gene + '.CCF', s]) or GSM.loc[gene + '.CCF', s] == -1:
            GSM.loc[gene + '.CCF', s] = event[ccf_col]
            if args.alt_count_fn:
                alt_count_df.loc[gene, s] = event['t_alt_count']
        elif not np.isnan(event[ccf_col]):
            GSM.loc[gene + '.CCF', s] = max(GSM.loc[gene + '.CCF', s], event[ccf_col])
            if event[ccf_col] > GSM.loc[gene + '.CCF', s]:
                alt_count_df.loc[gene, s] = event['t_alt_count']


GSM = GSM.round(4)
GSM = GSM.replace(-1, 0.0)

if args.ploidy:
    ploidy_df = pd.read_csv(args.ploidy, sep='\t', header=None, index_col=0)
    GSM.loc['PLOIDY'] = ploidy_df.loc[GSM.columns, 1].values

if args.purity:
    purity_df = pd.read_csv(args.purity, sep='\t', header=None, index_col=0)
    GSM.loc['PURITY'] = purity_df.loc[GSM.columns, 1].values

if args.coo:
    coo_df = pd.read_csv(args.coo, sep='\t', index_col=0)
    GSM.loc['COO'] = coo_df.loc[GSM.columns, 'COO'].values

if args.alt_count_fn:
    alt_count_df.index = alt_count_df.index.str.upper()
    alt_count_df.to_csv(args.alt_count_fn, sep='\t')

GSM.index = GSM.index.str.upper()
GSM.to_csv(args.output_fn, sep='\t')
