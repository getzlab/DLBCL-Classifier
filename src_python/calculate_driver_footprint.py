import pandas as pd
import re
import math

qval_cutoff = 0.10
datafile = "../data_tables/gsm/DLBCL.699.fullGSM.Sep_23_2022.tsv"
sv_bedfile = "../data_tables/feature_bp_counts/DLBCL_Rearrangm_probes_v1_AllTracks-REFORMAT.bed"
scnas_file = "../data_tables/feature_bp_counts/SCNAs_regions.tsv"
qval_file = '../data_tables/qval_dfs/fisher_exact_5x2.Sep_23_2022.combined.tsv'
arm_file = '../data_tables/feature_bp_counts/hg19.GISTIC.arms.tsv'
sig_genes = '../data_tables/mutsig2cv_gistic_qvalues/DLBCL_550_training_noPDE4DIP_noHISTartifacts.sig_genes.txt'

qvalDF = pd.read_csv(qval_file, sep='\t', index_col=0)
armDF = pd.read_csv(arm_file, sep='\t')
data = pd.read_csv(datafile, sep='\t', index_col=0)
svs = pd.read_csv(sv_bedfile, header=None, sep='\t')
scnas = pd.read_csv(scnas_file, sep='\t')
gene_regions = pd.read_csv(sig_genes, header=0, sep='\t', low_memory=False, index_col=0)
qvalDF = qvalDF.loc[qvalDF['q'] <= qval_cutoff]

armDF['arm'] = armDF['arm'].str.upper()
armDF = armDF.set_index(['arm'])

data = data.loc[qvalDF.index]

footprintSizes = {'IG_SVs': []}
full_bp_lengths = {'IG_SVs': []}
for gene in data.index:
    footprintSizes[gene] = []
    full_bp_lengths[gene] = []

########################
# SV & SCNA formatting #
########################

sv_filter = svs.iloc[:, 3].str.contains('BCL6') | svs.iloc[:, 3].str.contains('BCL2') | svs.iloc[:, 3].str.contains('MYC') | svs.iloc[:, 3].str.contains('IG')
svs = svs.loc[sv_filter, :]

scna_filter = data.index.str.contains('.AMP|.DEL')
scna_events_data = data.loc[scna_filter, :].index
scna_events_data = set(scna_events_data)

scnas['Name'] = 'X' + scnas['Name'].str.replace(':', '.')
scna_events_NatMed = scnas['Name']
scna_events_NatMed = set(scna_events_NatMed)

#############################
# Mutations footprint calcs #
#############################
for gene in data.index:
    gene2 = gene.replace('.', '-')
    if gene2 in gene_regions.index:
        footprintSizes[gene].append(gene_regions.loc[gene2, 'codelen'])
        full_bp_lengths[gene].append(gene_regions.loc[gene2, 'codelen'])

######################
# SV footprint calcs #
######################

for idx, row in svs.iterrows():
    if 'BCL6' in row[3]:
        val = 120 * math.ceil(float(row[2] - row[1] + 1) / 120)
        footprintSizes['SV.BCL6'].append(val)
        full_bp_lengths['SV.BCL6'].append(val)
    elif 'BCL2' in row[3]:
        val = 120 * math.ceil(float(row[2] - row[1] + 1) / 120)
        footprintSizes['SV.BCL2'].append(val)
        full_bp_lengths['SV.BCL2'].append(row[2] - row[1] + 1)
    elif 'MYC' in row[3]:
        val = 120 * math.ceil(float(row[2] - row[1] + 1) / 120)
        footprintSizes['SV.MYC'].append(val)
        full_bp_lengths['SV.MYC'].append(val)
    elif 'IG' in row[3]:
        val = 120 * math.ceil(float(row[2] - row[1] + 1) / 120)
        footprintSizes['IG_SVs'].append(val)
        full_bp_lengths['IG_SVs'].append(val)

########################
# SCNA footprint calcs #
########################

boundaries_reg = re.compile('^.+:(\d+)-(\d+)\(.*')
for idx, row in scnas.iterrows():
    curr_region = row['Name']
    if curr_region == 'X19Q13.32.DEL':
        curr_region = 'X19Q13.32.1.DEL'
    if curr_region == 'X18Q21.33(BCL2).AMP':
        curr_region = 'X18Q22.2.AMP'

    if curr_region not in footprintSizes:
            continue
    if row['Type'] == 'Arm' or row['Type'] == 'arm':
        armstart = armDF.loc[curr_region[1::].replace('.AMP', '').replace('.DEL', '')]['start']
        armend = armDF.loc[curr_region[1::].replace('.AMP', '').replace('.DEL', '')]['xEnd']
        numprobes = math.ceil(float(armend - armstart + 1) / 600000)
        armsize = numprobes * 120
        footprintSizes[curr_region].append(armsize)
        full_bp_lengths[curr_region].append(armend - armstart + 1)
    else:
        boundaries = boundaries_reg.findall(row['WidePeakLimits'])[0]
        numprobes = max(math.ceil(float(int(boundaries[1]) - int(boundaries[0]) + 1) / 200000), 10)
        focalsize = 120 * numprobes
        footprintSizes[curr_region].append(focalsize)
        full_bp_lengths[curr_region].append(int(boundaries[1]) - int(boundaries[0]) + 1)

gene_position = 0
targets_position = 0

# Need to map exomes to actual genes, so have to iterate through each dataframe and calculate
# if an exome region overlaps with a gene region.
# while targets_position < len(exome_targets.index) and gene_position < len(gene_regions.index):
#     exome_chr = exome_targets.iloc[targets_position]['chr']
#     exome_start = exome_targets.iloc[targets_position]['start']
#     exome_end = exome_targets.iloc[targets_position]['end']
#
#     gene_chr = gene_regions.iloc[gene_position]['#chrom']
#     gene_start = gene_regions.iloc[gene_position]['cdsStart']
#     gene_end = gene_regions.iloc[gene_position]['cdsEnd']
#
#     if gene_chr > exome_chr:
#         targets_position += 1
#         continue
#     if exome_chr > gene_chr:
#         gene_position += 1
#         continue
#     if ((exome_start <= gene_end) and (exome_start >= gene_start)) or \
#        ((exome_end <= gene_end) and (exome_end >= gene_start)) or \
#        ((gene_start <= exome_end) and (gene_start >= exome_start)) or \
#        ((gene_end <= exome_end) and (gene_end >= exome_start)):
#         if gene_regions.iloc[gene_position]['name2'] in data.index:
#             print(exome_chr, gene_chr, exome_start, gene_start, gene_start, gene_regions.iloc[gene_position]['name2'])
#             currSize = 120 * math.ceil(float(exome_end - exome_start + 1)/120)
#             footprintSizes[gene_regions.iloc[gene_position]['name2']].append(currSize)
#         targets_position += 1
#     elif exome_start >= gene_end:
#         gene_position += 1
#     elif gene_start >= exome_end:
#         targets_position += 1

for key, val in footprintSizes.items():
    footprintSizes[key] = sum(val)
    full_bp_lengths[key] = sum(full_bp_lengths[key])

# fix merged peak
footprintSizes['X18Q21.32.AMP'] = footprintSizes['X18Q22.2.AMP']

footprintDF = pd.DataFrame().from_dict(footprintSizes, orient='index')
full_bp_df = pd.DataFrame().from_dict(full_bp_lengths, orient='index')

IG_SV = footprintDF.loc['IG_SVs'][0]
BCL6_SV = footprintDF.loc['SV.BCL6'][0]
BCL2_SV = footprintDF.loc['SV.BCL2'][0]
MYC_SV = footprintDF.loc['SV.MYC'][0]

footprintDF.loc['IG_SVs'][0] = 0
footprintDF.loc['SV.BCL6'][0] = 0
footprintDF.loc['SV.BCL2'][0] = 0
footprintDF.loc['SV.MYC'][0] = 0
footprintDF = footprintDF.sort_values(by=[0])
footprintDF['integral'] = footprintDF[0].cumsum()
footprintDF['Total BP Count'] = full_bp_df[0]
footprintDF.loc['MYD88.L265P', 'Total BP Count'] = 1
footprintDF.loc['MYD88.OTHER', 'Total BP Count'] = footprintDF.loc['MYD88', 'Total BP Count'] - 1
footprintDF.columns = ['Classifier Footprint', 'integral', 'Total BP Count']
footprintDF.loc['MYD88.OTHER', 'integral'] = footprintDF.loc['MYD88', 'integral']
footprintDF.loc['MYD88.L265P', 'integral'] = footprintDF.loc['MYD88', 'integral']
footprintDF.index.name = 'Driver'
footprintDF.to_csv('../data_tables/feature_bp_counts/footprint_table_Sep_23_2022.tsv', sep='\t', header=True, index=True)

