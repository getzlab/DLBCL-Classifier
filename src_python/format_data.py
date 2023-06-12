import pandas as pd
from sklearn.decomposition import PCA
import numpy as np

# The target file must contain the columns named "C#" where # refers to the cluster number, and a column "cluster".
# So if 4 clusters, you should have 4 columns named C1, C2, C3, and C4.
# The value in each column, C#, refer to the sample's probability of C# membership

# Returns two pandas dataframes - a Sample x Feature input data frame, and a Sample x Cluster target data frame


def sigmoid(x):
    return 1/(1+np.exp(-1*x))


def weightfn(x):
    return 2*sigmoid(1/2*(x-10))+1


def construct_reduced_winning_version(data, add_missing_features=False):
    # This function expects samples in rows, features in columns
    if 'MYD88' in data.index:
        data = data.T

    qval_df = pd.read_csv('../data_tables/qval_dfs/fisher_exact_5x2.Sep_23_2022.combined.tsv', sep='\t', low_memory=False, index_col=0)
    genes_to_drop = list(data.columns[~data.columns.isin(qval_df.index)])
    data = data.drop(genes_to_drop, axis=1)

    if add_missing_features:
        missing_f = [x for x in qval_df.index if (x not in data.columns) and (qval_df.loc[x, 'q'] <= 0.10)]
        for feature in missing_f:
            print('Feature', feature, 'not found. Adding with full zeros.')
            data[feature] = 0

    data = data.astype(float).astype(int)

    # largest 5 removed manually from the vectors below
    # 'X1Q.AMP', 'X5Q.AMP', 'X4Q35.1.DEL', 'X1Q23.3.AMP', 'X9Q21.13.DEL'
    # C1 ##############

    list_bcl6 = ["SV.BCL6", "BCL6"]
    BCL6_ALT = data[list_bcl6].sum(axis=1)

    list_notch2_vec = ["NOTCH2", "SPEN", "DTX1"]
    NOTCH2_vec = data[list_notch2_vec].sum(axis=1)

    list_M88O_vec = ["MYD88.OTHER", "TNFAIP3", "TNIP1", "BCL10", "NFKBIE"]
    M88O_vec = data[list_M88O_vec].sum(axis=1)

    list_C1_vec4 = ["UBE2A", "TMEM30A", "ZEB2", "GNAI2", "X5P.AMP",
                    "POU2F2", "IKZF3", "X3Q28.DEL", "EBF1", "LYN", "BCL7A", "CXCR4",
                    "CCDC27", "TUBGCP5", "SMG7", "RHOA", "BTG2"]
    C1_vec4 = data[list_C1_vec4].sum(axis=1)

    list_CD70_vec = ["CD70", "FAS", "CD58", "B2M", "FADD", "HLA.B"]
    CD70_vec = data[list_CD70_vec].sum(axis=1)

    # C2 ##############

    list_tp53 = ["TP53", "X17P.DEL"]
    TP53_biallelic = data[list_tp53].sum(axis=1)

    X21Q_AMP = data["X21Q.AMP"]

    list_C2_sum_arm = ["X17P.DEL", "X21Q.AMP", "X11Q.AMP", "X6P.AMP", "X11P.AMP",
                       "X6Q.DEL", "X7P.AMP", "X13Q.AMP", "X7Q.AMP", "X3Q.AMP", "X5P.AMP",
                       "X18P.AMP", "X3P.AMP", "X19Q.AMP", "X9Q.AMP",
                       "X12P.AMP", "X12Q.AMP"]
    Sum_C2_ARM = data[list_C2_sum_arm].sum(axis=1)

    list_C2_sum_focal = ["X1P36.11.DEL", "X1P31.1.DEL", "X1P13.1.DEL", "X2Q22.2.DEL", "X16Q12.1.DEL",
                         "X14Q32.31.DEL", "X1P36.32.DEL", "X15Q15.3.DEL",
                         "X4Q21.22.DEL", "X9P21.3.DEL", "X8Q24.22.AMP", "X12P13.2.DEL", "X2P16.1.AMP",
                         "X8Q12.1.DEL", "X19P13.2.DEL", "X17Q25.1.DEL", "X1Q42.12.DEL", "X3P21.31.DEL",
                         "X18Q23.DEL", "X19P13.3.DEL", "X13Q34.DEL", "X7Q22.1.AMP", "X10Q23.31.DEL",
                         "X9P24.1.AMP", "X3Q28.AMP", "X11Q23.3.AMP", "X17Q24.3.AMP", "X3Q28.DEL",
                         "X13Q14.2.DEL", "X18Q21.32.AMP", "X19Q13.32.DEL", "X6P21.1.AMP", "X18Q22.2.AMP",
                         "EP300", "ZNF423", "CD274"]
    Sum_C2_FOCAL = data[list_C2_sum_focal].sum(axis=1)

    # C3 ##############

    list_bcl2 = ["BCL2", "SV.BCL2"]
    BCL2_combined = data[list_bcl2].sum(axis=1)

    list_CREBBP = ["CREBBP", "EZH2", "KMT2D", "EP300"]
    CREBBP_vec = data[list_CREBBP].sum(axis=1)

    list_GNA13 = ["GNA13", "TNFRSF14", "MAP2K1", "MEF2B", "IRF8", "HVCN1",
                  "GNAI2", "MEF2C", "SOCS1", "EEF1A1", "RAC2", "X12Q.AMP",
                  "POU2AF1", "X6Q14.1.DEL"]
    GNA13_vec = data[list_GNA13].sum(axis=1)

    PTEN = data["PTEN"] + data["X10Q23.31.DEL"] + data["X13Q14.2.DEL"]

    SV_MYC = data["SV.MYC"]

    # C4 ##############

    list_Hist_comp = ["HIST1H2AC", "HIST1H1E", "HIST1H1B", "HIST1H2AM", "HIST1H1C", "HIST1H1D", "HIST1H2BC"]
    Hist_comp = data[list_Hist_comp].sum(axis=1)

    list_SGK1_vec = ["SGK1", "TET2", "NFKBIA", "STAT3", "PTPN6", "BRAF", "KRAS",
                     "CD83", "SF3B1", "CD274", "MEF2C", "KLHL6", "CXCR4", "PTEN",
                     "RAC2", "SESN3", "SOCS1", "METAP1D"]
    SGK1_vec = data[list_SGK1_vec].sum(axis=1)

    list_DUSP2_vec = ["DUSP2", "ZFP36L1", "CRIP1", "ACTB", "LTB", "YY1", "PABPC1"]
    DUSP2_vec = data[list_DUSP2_vec].sum(axis=1)

    # C5 ##############

    list_TBL1XR1_vec = ["TBL1XR1", "PIM1", "PRDM1", "ETV6", "ZC3H12A", "BTG1", "BTG2",
                        "IGLL5", "TMSB4X", "GRHPR", "HLA.C", "MYD88", "TOX", "LYN",
                        "POU2F2", "IKZF3", "HLA.A", "ZFP36L1", "CARD11", "SF3B1",
                        "HLA.B", "IRF2BP2", "OSBPL10", "ATP2A2", "PIM2", "IRF4", "BCL11A",
                        "METAP1D", "ETS1", "CCDC27"]
    TBL1XR1_vec = data[list_TBL1XR1_vec].sum(axis=1)

    MYD88_L265P_CD79B = data["MYD88.L265P"] + data["CD79B"]

    list_Sum_C5_sig = ["X18Q.AMP", "X3Q.AMP", "X3P.AMP", "X19Q13.42.AMP", "X6Q21.DEL",
                       "X18P.AMP", "X19Q.AMP", "X8Q12.1.DEL", "X6Q14.1.DEL", "X19P13.2.DEL",
                       "X9P21.3.DEL", "X18Q21.32.AMP", "X18Q22.2.AMP", "X1Q42.12.DEL", "X1Q32.1.AMP", "X6P21.33.DEL"]
    Sum_C5_CNA = data[list_Sum_C5_sig].sum(axis=1)

    # MISC ############
    CN_2P16_1_AMP = data["X2P16.1.AMP"]

    reduced_feature_dict = \
        {'BCL6_ALT': BCL6_ALT,
         'NOTCH2_vec': NOTCH2_vec,
         'M88O_vec': M88O_vec,
         'C1_vec4': C1_vec4,
         'CD70_vec': CD70_vec,
         'TP53_biallelic': TP53_biallelic,
         'X21Q_AMP': X21Q_AMP,
         'Sum_C2_ARM': Sum_C2_ARM,
         'Sum_C2_FOCAL': Sum_C2_FOCAL,
         'BCL2_combined': BCL2_combined,
         'CREBBP_vec': CREBBP_vec,
         'GNA13_vec': GNA13_vec,
         'PTEN': PTEN,
         'SV_MYC': SV_MYC,
         'Hist_comp': Hist_comp,
         'SGK1_vec': SGK1_vec,
         'DUSP2_vec': DUSP2_vec,
         'CN_2P16_1_AMP': CN_2P16_1_AMP,
         'TBL1XR1_vec': TBL1XR1_vec,
         'MYD88_L265P_CD79B': MYD88_L265P_CD79B,
         'Sum_C5_CNA': Sum_C5_CNA}

    data = pd.DataFrame.from_dict(reduced_feature_dict)

    return data


def format_inputs(datafile, target_file, training_set,
                  no_sv=False,
                  no_cna=False,
                  ploidy=False,
                  no_silent=False,
                  use_pca=False,
                  coo=False,
                  reduced_version=None,
                  no_muts=False,
                  qval=0.10,
                  no_arms=False,
                  no_focals=False,
                  bp_cutoff=None,
                  drop_empty_vectors=True,
                  remove_outliers=False,
                  weightl265p=False,
                  remove_largest_n=None,
                  nosvbcl6=False):

    if reduced_version and reduced_version != '3.4':
        print('Wrong version, should be reduced 3.4. Exiting.')
        exit()

    data = pd.read_csv(datafile, delimiter='\t', index_col=0, low_memory=False)
    qval_df = pd.read_csv('../data_tables/qval_dfs/fisher_exact_5x2.Sep_23_2022.combined.tsv', sep='\t', low_memory=False, index_col=0)
    targets = pd.read_csv(target_file, delimiter='\t', index_col=0)

    if 'Driver_Discovery' in data.columns:
        data.drop(['Driver_Discovery'], axis=1, inplace=True)

    data = data[training_set]
    data = data.transpose()

    if weightl265p:
        data['MYD88.L265P'] = data['MYD88.L265P'].replace('2.0', '3.0')

    # Set genome doubling & drop ploidy
    # Added in at the end if ploidy is flagged
    if 'PLOIDY' in data.columns:
        genome_doubling = data['PLOIDY'].astype(float).copy(deep=True)
        genome_doubling[genome_doubling >= 3.00] = 3
        genome_doubling[genome_doubling < 3.00] = 0
        data.drop('PLOIDY', axis=1, inplace=True)

    if ploidy:
        data['genome_doubling'] = genome_doubling

    # COO block, numerical distance-based encoding (vector of 3, ABC/GCB/UNC). NA samples must be dropped.
    # These will be added in at the end if coo is flagged
    if coo:
        data = data.loc[data['COO'] != 'na']
        coo_col = data['COO'].map({'ABC': 2, 'GCB': 0, 'UNC': 1})

    # Zero out the largest N features
    if remove_largest_n is not None:
        bptable = pd.read_csv('../data_tables/feature_bp_counts/footprint_table_Sep_23_2022.tsv', sep='\t', index_col=0)
        bptable.index = bptable.index.str.replace('-', '.')
        features_to_remove = bptable.tail(remove_largest_n).index
        features_to_remove_c = []
        for fr in features_to_remove:
            if 'SV' in fr:
                continue
            if fr in data.columns:
                data[fr] = 0
                features_to_remove_c.append(fr)
            elif 'X'+fr+'.AMP' in data.columns:
                data['X'+fr+'.AMP'] = 0
                features_to_remove_c.append('X'+fr+'.AMP')
            elif 'X'+fr+'.DEL' in data.columns:
                data['X'+fr+'.DEL'] = 0
                features_to_remove_c.append('X'+fr+'.DEL')
            else:
                print('Feature', fr,'not found, exiting!')
                exit()

        fn = '../data_tables/feature_bp_counts/trim_n' + str(remove_largest_n) + '_removedfeatures.tsv'
        with open(fn, 'w+') as f:
            for feature in features_to_remove_c:
                f.write(feature + '\n')

        keptfeatures = sorted(list(set(bptable.index) - set(features_to_remove_c)))
        fn = '../data_tables/feature_bp_counts/trim_n' + str(remove_largest_n) + '_keptfeatures.tsv'
        with open(fn, 'w+') as f:
            for feature in keptfeatures:
                f.write(feature + '\n')

    # Zero out features that exceed footprint quota
    if bp_cutoff is not None:
        bptable = pd.read_csv('../data_tables/feature_bp_counts/footprint_table_Sep_23_2022.tsv', sep='\t', index_col=0)
        bptable.index = bptable.index.str.replace('-', '.')
        features_to_remove = bptable['integral'][bptable['integral'] >= bp_cutoff].index
        features_to_remove = [x.replace('-', '.') for x in features_to_remove]
        features_to_remove_c = []
        for fr in features_to_remove:
            if 'SV' in fr:
                continue
            count = 0
            if fr in data.columns:
                data[fr] = 0
                features_to_remove_c.append(fr)
                count += 1
            if 'X'+fr+'.AMP' in data.columns:
                data['X'+fr+'.AMP'] = 0
                features_to_remove_c.append('X'+fr+'.AMP')
                count += 1
            if 'X'+fr+'.DEL' in data.columns:
                data['X'+fr+'.DEL'] = 0
                features_to_remove_c.append('X'+fr+'.DEL')
                count += 1
            if not count:
                print('Feature',  fr, 'not found! Exiting!')
                exit()

        fn = '../data_tables/feature_bp_counts/cutoff_' + str(bp_cutoff) + '_removedfeatures.tsv'
        with open(fn, 'w+') as f:
            for feature in features_to_remove_c:
                f.write(feature+'\n')

        keptfeatures = sorted(list(set(bptable.index) - set(features_to_remove_c)))
        fn = '../data_tables/feature_bp_counts/cutoff_' + str(bp_cutoff) + '_keptfeatures.tsv'
        with open(fn, 'w+') as f:
            for feature in keptfeatures:
                f.write(feature+'\n')

    # Zero out other features that we should drop
    if no_sv:
        sv_cols = data.columns[data.columns.str.contains('SV.', regex=False)]
        data[sv_cols] = 0

    if nosvbcl6:
        data['SV.BCL6'] = 0

    if no_muts:
        sv_cols = data.columns[data.columns.str.contains('SV.', regex=False)]
        cna_cols = data.columns[data.columns.str.contains('.AMP', regex=False) | data.columns.str.contains('.DEL', regex=False)]
        mut_cols = data.columns[~(data.columns.isin(sv_cols) | data.columns.isin(cna_cols))]
        data[mut_cols] = 0

    # no focals block - set all focal cn events to 0
    if no_focals or no_cna:
        allarms = []
        for i in range(1, 23):
            allarms.append('X' + str(i) + 'P.DEL')
            allarms.append('X' + str(i) + 'Q.DEL')
            allarms.append('X' + str(i) + 'P.AMP')
            allarms.append('X' + str(i) + 'Q.AMP')
        for feature in data.columns:
            if feature.endswith('.DEL') | feature.endswith('.AMP'):
                if feature not in allarms:
                    data[feature] = 0

    # no arms block - set all arm events to 0
    if no_arms or no_cna:
        allarms = []
        for i in range(1, 23):
            allarms.append('X' + str(i) + 'P.DEL')
            allarms.append('X' + str(i) + 'Q.DEL')
            allarms.append('X' + str(i) + 'P.AMP')
            allarms.append('X' + str(i) + 'Q.AMP')
        for feature in data.columns:
            if feature in allarms:
                data[feature] = 0

    # only drop all features if we aren't using reduced recipe, since reduced will expect all features
    if not reduced_version:
        data = data[qval_df.loc[qval_df['q'] <= qval].index]
    else:
        data = data[qval_df.loc[qval_df['q'] <= 0.10].index]

    targets = targets.loc[targets.index.isin(data.index)]
    data = data.loc[targets.index]
    data = data.astype(float)

    if use_pca and not reduced_version:
        if no_sv:
            data = data[data.columns.drop(list(data.filter(regex='SV\.')))]
        if no_cna:
            data = data[data.columns.drop(list(data.filter(regex='\.DEL')))]
            data = data[data.columns.drop(list(data.filter(regex='\.AMP')))]
        pca = PCA(n_components=use_pca)
        p_components = pca.fit_transform(data)
        p_components = pd.DataFrame(p_components)
        p_components.index = data.index
        data = p_components

    if reduced_version == '3.4':
        if qval <= 0.10:
            genes_to_drop = list(qval_df.loc[qval_df['q'] > qval].index)
            genes_to_drop = [x for x in genes_to_drop if x in data.columns]
            genes_to_keep = [x for x in data.columns if x not in genes_to_drop]
            data[genes_to_drop] = 0

            fn = '../data_tables/feature_bp_counts/reduced_cutoff_qval' + str(qval) + '_removedfeatures.tsv'
            with open(fn, 'w+') as f:
                f.write('\n'.join(genes_to_drop))

            fn = '../data_tables/feature_bp_counts/reduced_cutoff_qval' + str(qval) + '_keptfeatures.tsv'
            with open(fn, 'w+') as f:
                f.write('\n'.join(genes_to_keep))

        # C1 ##############

        list_bcl6 = ["SV.BCL6", "BCL6"]
        BCL6_ALT = data[list_bcl6].sum(axis=1)

        list_notch2_vec = ["NOTCH2", "SPEN", "DTX1"]
        NOTCH2_vec = data[list_notch2_vec].sum(axis=1)

        list_M88O_vec = ["MYD88.OTHER", "TNFAIP3", "TNIP1", "BCL10", "NFKBIE"]
        M88O_vec = data[list_M88O_vec].sum(axis=1)

        list_C1_vec4 = ["UBE2A", "TMEM30A", "ZEB2", "GNAI2", "X5P.AMP", "X5Q.AMP",
                        "POU2F2", "IKZF3", "X3Q28.DEL", "EBF1", "LYN", "BCL7A", "CXCR4",
                        "CCDC27", "TUBGCP5", "SMG7", "RHOA", "BTG2"]
        C1_vec4 = data[list_C1_vec4].sum(axis=1)

        list_CD70_vec = ["CD70", "FAS", "CD58", "B2M", "FADD", "HLA.B"]
        CD70_vec = data[list_CD70_vec].sum(axis=1)

        # C2 ##############

        list_tp53 = ["TP53", "X17P.DEL"]
        TP53_biallelic = data[list_tp53].sum(axis=1)

        X21Q_AMP = data["X21Q.AMP"]

        list_C2_sum_arm = ["X17P.DEL", "X21Q.AMP", "X11Q.AMP", "X6P.AMP", "X11P.AMP",
                           "X6Q.DEL", "X7P.AMP", "X13Q.AMP", "X7Q.AMP", "X3Q.AMP", "X5P.AMP",
                           "X5Q.AMP", "X18P.AMP", "X3P.AMP", "X19Q.AMP", "X9Q.AMP", "X1Q.AMP",
                           "X12P.AMP", "X12Q.AMP"]
        Sum_C2_ARM = data[list_C2_sum_arm].sum(axis=1)

        list_C2_sum_focal = ["X1P36.11.DEL", "X1P31.1.DEL", "X1P13.1.DEL", "X2Q22.2.DEL", "X16Q12.1.DEL",
                             "X14Q32.31.DEL", "X1P36.32.DEL", "X4Q35.1.DEL", "X9Q21.13.DEL", "X15Q15.3.DEL",
                             "X4Q21.22.DEL", "X9P21.3.DEL", "X8Q24.22.AMP", "X12P13.2.DEL", "X2P16.1.AMP",
                             "X8Q12.1.DEL", "X19P13.2.DEL", "X17Q25.1.DEL", "X1Q42.12.DEL", "X3P21.31.DEL",
                             "X18Q23.DEL", "X19P13.3.DEL", "X13Q34.DEL", "X7Q22.1.AMP", "X10Q23.31.DEL",
                             "X9P24.1.AMP", "X1Q23.3.AMP", "X3Q28.AMP", "X11Q23.3.AMP", "X17Q24.3.AMP", "X3Q28.DEL",
                             "X13Q14.2.DEL", "X18Q21.32.AMP", "X19Q13.32.DEL", "X6P21.1.AMP", "X18Q22.2.AMP",
                             "EP300", "ZNF423", "CD274"]
        Sum_C2_FOCAL = data[list_C2_sum_focal].sum(axis=1)

        # C3 ##############

        list_bcl2 = ["BCL2", "SV.BCL2"]
        BCL2_combined = data[list_bcl2].sum(axis=1)

        list_CREBBP = ["CREBBP", "EZH2", "KMT2D", "EP300"]
        CREBBP_vec = data[list_CREBBP].sum(axis=1)

        list_GNA13 = ["GNA13", "TNFRSF14", "MAP2K1", "MEF2B", "IRF8", "HVCN1",
                      "GNAI2", "MEF2C", "SOCS1", "EEF1A1", "RAC2", "X12Q.AMP",
                      "POU2AF1", "X6Q14.1.DEL"]
        GNA13_vec = data[list_GNA13].sum(axis=1)

        PTEN = data["PTEN"] + data["X10Q23.31.DEL"] + data["X13Q14.2.DEL"]

        SV_MYC = data["SV.MYC"]

        # C4 ##############

        list_Hist_comp = ["HIST1H2AC", "HIST1H1E", "HIST1H1B", "HIST1H2AM", "HIST1H1C", "HIST1H1D", "HIST1H2BC"]
        Hist_comp = data[list_Hist_comp].sum(axis=1)

        list_SGK1_vec = ["SGK1", "TET2", "NFKBIA", "STAT3", "PTPN6", "BRAF", "KRAS",
                         "CD83", "SF3B1", "CD274", "MEF2C", "KLHL6", "CXCR4", "PTEN",
                         "RAC2", "SESN3", "SOCS1", "METAP1D"]
        SGK1_vec = data[list_SGK1_vec].sum(axis=1)

        list_DUSP2_vec = ["DUSP2", "ZFP36L1", "CRIP1", "ACTB", "LTB", "YY1", "PABPC1"]
        DUSP2_vec = data[list_DUSP2_vec].sum(axis=1)

        # C5 ##############

        list_TBL1XR1_vec = ["TBL1XR1", "PIM1", "PRDM1", "ETV6", "ZC3H12A", "BTG1", "BTG2",
                            "IGLL5", "TMSB4X", "GRHPR", "HLA.C", "MYD88", "TOX", "LYN",
                            "POU2F2", "IKZF3", "HLA.A", "ZFP36L1", "CARD11", "SF3B1",
                            "HLA.B", "IRF2BP2", "OSBPL10", "ATP2A2", "PIM2", "IRF4", "BCL11A",
                            "METAP1D", "ETS1", "CCDC27"]
        TBL1XR1_vec = data[list_TBL1XR1_vec].sum(axis=1)

        MYD88_L265P_CD79B = data["MYD88.L265P"] + data["CD79B"]

        list_Sum_C5_sig = ["X18Q.AMP", "X3Q.AMP", "X3P.AMP", "X19Q13.42.AMP", "X6Q21.DEL",
                           "X18P.AMP", "X19Q.AMP", "X8Q12.1.DEL", "X6Q14.1.DEL", "X19P13.2.DEL",
                           "X9P21.3.DEL", "X18Q21.32.AMP", "X18Q22.2.AMP", "X1Q42.12.DEL", "X1Q32.1.AMP", "X6P21.33.DEL"]
        Sum_C5_CNA = data[list_Sum_C5_sig].sum(axis=1)

        # MISC ############
        CN_2P16_1_AMP = data["X2P16.1.AMP"]

        reduced_feature_dict = \
            {'BCL6_ALT': BCL6_ALT,
             'NOTCH2_vec': NOTCH2_vec,
             'M88O_vec': M88O_vec,
             'C1_vec4': C1_vec4,
             'CD70_vec': CD70_vec,
             'TP53_biallelic': TP53_biallelic,
             'X21Q_AMP': X21Q_AMP,
             'Sum_C2_ARM': Sum_C2_ARM,
             'Sum_C2_FOCAL': Sum_C2_FOCAL,
             'BCL2_combined': BCL2_combined,
             'CREBBP_vec': CREBBP_vec,
             'GNA13_vec': GNA13_vec,
             'PTEN': PTEN,
             'SV_MYC': SV_MYC,
             'Hist_comp': Hist_comp,
             'SGK1_vec': SGK1_vec,
             'DUSP2_vec': DUSP2_vec,
             'CN_2P16_1_AMP': CN_2P16_1_AMP,
             'TBL1XR1_vec': TBL1XR1_vec,
             'MYD88_L265P_CD79B': MYD88_L265P_CD79B,
             'Sum_C5_CNA': Sum_C5_CNA}

        data = pd.DataFrame.from_dict(reduced_feature_dict)

    # Add in GENOME_DOUBLING and/or COO based on flags
    if ploidy:
        genome_doubling = genome_doubling[data.index]
        data['GENOME_DOUBLING'] = genome_doubling

    if coo:
        coo_col = coo_col[data.index]
        data['COO'] = coo_col

    if drop_empty_vectors:
        data = data[(data.sum(axis=0) != 0)[(data.sum(axis=0) != 0)].index]

    return data, targets[['C1', 'C2', 'C3', 'C4', 'C5']]
