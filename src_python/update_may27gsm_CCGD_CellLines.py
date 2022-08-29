import pandas as pd


# map the cell lines

mappedCLs = {'Balm3_L_156272': 'BALM3',
             'BL-2_L_156270': 'BL_2',
             'BL-41_L_156263': 'BL_41',
             'CA-46_L_156268': 'CA_46',
             'CTB-1_L_156261': 'CTB_1',
             'DB_L_156242': 'DB',
             'DHL2': 'DHL2',
             'DHL10_L_156238': 'DHL10',
             'DHL16_L_156271': 'DHL16',
             'DHL4_L_156223': 'DHL4',
             'DHL5_L_156273': 'DHL5',
             'DHL6_L_156224': 'DHL6',
             'DHL7_L_156236': 'DHL7',
             'DHL8_L_156237': 'DHL8',
             'Daudi_L_156265': 'DAUDI',
             'DoHH2_L_156259': 'DOHH2',
             'Dogkit_L_156269': 'DOGKIT',
             'FARAGE_L_156248': 'FARAGE',
             'GUMBUS_L_156256': 'GUMBUS',
             'HBL1_L_156227': 'HBL1',
             'HDLM2_L_156252': 'HDLM2',
             'HKBML_L_156257': 'HKBML',
             'HT_L_156243': 'HT',
             'K422_L_156230': 'K422',
             'KMH2_L_156251': 'KMH2',
             'Karpas1106K1106_L_156247': 'KARPAS1106',
             'L1236_L_156254': 'L1236',
             'L428_L_156250': 'L428',
             'L540_L_156253': 'L540',
             'Ly1_L_156225': 'LY1',
             'Ly10_L_156241': 'LY10',
             'Ly18_L_156239': 'LY18',
             'Ly19_L_156234': 'LY19',
             'LY3_L_156235': 'LY3',
             'Ly4_L_156231': 'LY4',
             'Ly7_L_156226': 'LY7',
             'Ly8_L_156240': 'LY8',
             'NU-DUL-1_L_156258': 'NU_DUL_1',
             'Namalwa_L_156266': 'NAMALWA',
             'Pfeiffer_L_156233': 'PFEIFFER',
             'Raji_L_156267': 'RAJI',
             'Ramos_L_156264':  'RAMOS',
             'SC1_L_156246': 'SC1',
             'SUPHDI_L_156249': 'SUPHDI',
             'TK_L_156255': 'TK',
             'TMD8_L_156228': 'TMD8',
             'TOLEDO_L_156232': 'TOLEDO',
             'U-H01_L_156262': 'U_H01',
             'U2932_L_156229': 'U2932',
             'U2940': 'U2940',
             'WSU-DLBCL2_L_156245': 'WSU_DLBCL2',
             'WSU-FSCCL_L_156260': 'WSU_FSCCL',
             'WSU-NHL_L_156244': 'WSU_NHL'}


datafile = '../data_tables/gsm/CellLines_addedgenesTET2.20-Apr-2020.tsv'
qval_file = '../data_tables/qval_dfs/fisher_exact_5x2_17-Aug-2022.combined.tsv'

cl_maf = '../data_tables/maf_files/CCGD_Lymphoma_Cell_Lines_53.with_noncoding.maf.annotated.trimmed'
cl_maf = pd.read_csv(cl_maf, sep='\t', low_memory=False)
#print('Loaded CL maf')

cell_lines = [x.upper() for x in mappedCLs.values()]

ubermatrix = pd.read_csv(datafile, sep='\t', index_col=0)
ubermatrix.index = ubermatrix.index.str.upper()
qval_df = pd.read_csv(qval_file, sep='\t', index_col=0)
qval_df = qval_df.loc[qval_df['q'] <= 0.10]


sv_indels = [x for x in ubermatrix.index if '.AMP' in x or '.DEL' in x or 'SV.' in x]
ubermatrix = ubermatrix.loc[sv_indels]

genes_to_add = [x for x in qval_df.index if x not in ubermatrix.index]
genes_to_add = [x.replace('.', '-') for x in genes_to_add]

cell_lines = [x.replace('-', '_') for x in cell_lines]
cell_lines_matrix = ubermatrix[cell_lines]

for gene in genes_to_add:
    appDF = pd.DataFrame([[0]*len(cell_lines_matrix.columns)], columns=cell_lines_matrix.columns, index=[gene])
    cell_lines_matrix = cell_lines_matrix.append(appDF)


noncoding_types = {'Intron', "5'UTR", "3'UTR", "5'Flank", "3'Flank", 'IGR', 'COULD_NOT_DETERMINE', 'RNA'}

for idx, row in cl_maf.iterrows():
    gene = row['Hugo_Symbol']
    if gene not in genes_to_add:
        continue

    currCL = row['Tumor_Sample_Barcode']
    vc = row['Variant_Classification']

    if row['Variant_Classification'] in noncoding_types:
        continue

    if currCL not in mappedCLs:
        continue

    mapped_cl = mappedCLs[currCL]

    if row['Variant_Classification'] == "5'Flank":
        continue

    if gene == 'MYD88':
        if row['Start_position'] == 38182641 and row['Tumor_Seq_Allele2'] == 'C':
            cell_lines_matrix.loc['MYD88-L265P', mapped_cl] = 2
        else:
            if vc == 'Silent':
                if cell_lines_matrix.loc['MYD88-OTHER', mapped_cl] == 0:
                    cell_lines_matrix.loc['MYD88-OTHER', mapped_cl] = 1
            else:
                cell_lines_matrix.loc['MYD88-OTHER', mapped_cl] = 2
    else:
        if vc == 'Silent':
            val = 1
        else:
            val = 2
        cell_lines_matrix.at[gene, mapped_cl] = max(val, cell_lines_matrix.loc[gene, mapped_cl])

cell_lines_matrix.index = [x.replace('-', '.') for x in cell_lines_matrix.index]
cell_lines_matrix.loc['MYD88'] = cell_lines_matrix.loc[['MYD88.L265P', 'MYD88.OTHER']].max(axis=0)
cell_lines_matrix.to_csv('../data_tables/gsm/GSM.CCGD.updated.Aug-17-2022.tsv', sep='\t')
with open('../data_tables/gsm/CCGD_Cell_lines.txt', 'w+') as f:
    for cl in cell_lines:
        f.write(cl+'\n')
