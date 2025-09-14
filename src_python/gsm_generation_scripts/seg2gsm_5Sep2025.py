import pandas as pd
import matlab_functions as mf
import numpy as np
import datetime
import argparse
import time
from datetime import datetime
TODAY = datetime.now().strftime("%d%b%Y")


parser = argparse.ArgumentParser()
parser.add_argument('-i', '--id',
                    help='Cohort set name.',
                    required=False, default="DLBCL.NCI_10", type=str)
parser.add_argument('-s', '--sample_set',
                    help='Sample file listing all samples\' in cohort, even w/o CNV event.',
                    required=False, default="", type=str)
parser.add_argument('-v', '--seg_file',
                    help='The seg file (CBS/IGV format) to find your samples\' CNV events.',
                    required=False, default="/private/tmp/CopyNumber_476_NCICCR_phs001444.log2CR_NCI_10_DLBCL10450-DLBCL10460.DNAcopy.1.70.0.seg",type=str)
#                   required=False, default="/private/tmp/DLBCL_699_segs.2021-12-15.NCI_10.06Sep2025.seg",type=str)
parser.add_argument('-x', '--cnv_blacklist_file',
                    help='CNV regions to exclude.',
                    required=False, type=str, default='../../data_tables/additional_gsm_inputs/CNV.hg19.bypos.111213.CR1_event_added.bed')
parser.add_argument('-b', '--broad_significance_file',
                    help='CNV arm-level regions to include.',
                    required=False, type=str, default='../../data_tables/additional_gsm_inputs/DLBCL_cHL_PMBL_broad_significance.txt')
parser.add_argument('-f', '--focal_file',
                    help='CNV foca; regions to include.',
                    required=False, type=str,default='../../data_tables/additional_gsm_inputs/all_significant_focal_peaks_with_cohort.txt')
parser.add_argument('-o', '--output_dir',
                    help='Output directory.',
                    required=False, type=str,default='/tmp/')
parser.add_argument('-g','--genome_build',
                    help='Genome build: hg19, hg38.',
                    required=False, type=str, default='hg19')
parser.add_argument('-d','--debug',
                    help='if not blank print intermediate output with debug string stub',
                    required=False, type=str, default='')


args = parser.parse_args()
seg_file = args.seg_file
cnv_blacklist_file = args.cnv_blacklist_file
broad_significance_file = args.broad_significance_file
focal_file = args.focal_file
output_dir = args.output_dir
debug = args.debug
OUTPUT_FN = args.output_dir + args.id + '.' + TODAY + '.gsm.tsv'

"""
OUTPUT_FN = '../../data_tables/gsm/DLBCL.699.scnaGSM.Sep_23_2022.tsv'
OUTPUT_FN = '../../data_tables/gsm/DLBCL.699.scnaGSM.Sep_23_2022.05Sep2025.tsv'

seg_file = '../../data_tables/additional_gsm_inputs/DLBCL_786_segs.2021-12-15.ccf.txt'

cnv_blacklist_file = '../../data_tables/additional_gsm_inputs/CNV.hg19.bypos.111213.CR1_event_added.bed'
broad_significance_file = '../../data_tables/additional_gsm_inputs/DLBCL_cHL_PMBL_broad_significance.txt'
focal_file = '../../data_tables/additional_gsm_inputs/all_significant_focal_peaks_with_cohort.txt'
"""
segs = pd.read_csv(seg_file, sep='\t')

single_amp_threshold = 0.1
single_del_threshold = -0.1
double_amp_threshold = 0.9
double_del_threshold = -0.9
arm_length_fraction_threshold = 2

AL = pd.read_csv(focal_file, sep='\t')
arm_level_significance = pd.read_csv(broad_significance_file, sep='\t')

# simple DNAcopy seg file format 
# ID	chrom	loc.start	loc.end	num.mark	seg.mean
# Sample Chromosome	Start	End	...	Segment_<Man
c = segs.columns
segs.rename(columns={c[0]: 'Sample',c[1]:'Chromosome',c[2]:'Start',c[3]:'End',c[-1]:'Segment_Mean'}, inplace=True)

# remove X and Y - seg2gsm file is for autosomes
segs = segs.loc[~segs['Chromosome'].isin(['X', 'Y'])].copy(deep=True)

# genomic position in 1D 
#segs.loc[:, 'gstart'] = mf.xhg19(segs['Chromosome'], segs['Start'])
#segs.loc[:, 'gend'] = mf.xhg19(segs['Chromosome'], segs['End'])

if args.genome_build == 'hg19':
    segs.loc[:, 'gstart'] = mf.xhg19(segs['Chromosome'], segs['Start'])
    segs.loc[:, 'gend'] = mf.xhg19(segs['Chromosome'], segs['End'])
if args.genome_build == 'hg38':
    segs.loc[:, 'gstart'] = mf.xhg38(segs['Chromosome'], segs['Start'])
    segs.loc[:, 'gend'] = mf.xhg38(segs['Chromosome'], segs['End'])


segs.loc[:, 'length'] = segs['End'] - segs['Start']

cnv_blacklist = pd.read_csv(cnv_blacklist_file, sep='\t')
segs = mf.apply_cnv_blacklist(segs, cnv_blacklist, AL, arm_level_significance)

# Segment_Mean in log2CR units 
segs['orignial_Segment_Mean'] = segs['Segment_Mean'].copy(deep=True)
for samp in sorted(segs['Sample'].unique()):
    sampseg = segs.loc[segs['Sample'] == samp].copy(deep=True)
    # sample_median in log2CR units 
    sample_median = mf.calc_region_median(sampseg, min(arm_level_significance['x1']), max(arm_level_significance['x2']), 2)
    segs.loc[segs['Sample'] == samp, 'sample_median'] = sample_median
    # recenter Segment_Mean (=log2CR) to median 
    segs.loc[segs['Sample'] == samp, 'Segment_Mean'] - segs.loc[segs['Sample'] == samp, 'Segment_Mean']  - sample_median

# save input Segment_Mean (=logCR) 
segs['log_segment_mean'] = segs['Segment_Mean'].copy(deep=True)
# convert Segment_Mean to 2*CR -2 units
segs['Segment_Mean'] = np.power(2, (segs['Segment_Mean'] + 1)) - 2.0



if len(debug)   > 0:
    f1 = output_dir + debug + '.segs.after_scaling.' + TODAY + '.seg'
    segs.to_csv(f1, sep='\t', index=False)

BA = arm_level_significance.loc[(arm_level_significance['significant_amplification'] == 1) &
                                (arm_level_significance['amplification_cohort'].str.contains('DLBCL'))]

BD = arm_level_significance.loc[(arm_level_significance['significant_deletion'] == 1) &
                                (arm_level_significance['deletion_cohort'].str.contains('DLBCL'))]

# Arm dels
arm_del_df = pd.DataFrame(0, index=BD['arm'] + '.DEL', columns=sorted(segs['Sample'].unique()))
#arm_del_df_ccf = arm_del_df.copy(deep=True)
#arm_del_df_ccf.index = arm_del_df_ccf.index + '.CCF'
#arm_del_df = pd.concat([arm_del_df, arm_del_df_ccf])

for patient in arm_del_df.columns:
    for arm in BD['arm']:

        boundix = (BD['arm'] == arm)
        boundix = boundix[boundix].index

        bd_x1 = BD.loc[boundix, 'x1'].values[0]
        bd_x2 = BD.loc[boundix, 'x2'].values[0]

        # (
        # (segs.gstart <= BD.x1(boundix) & segs.gend >= BD.x1(boundix)) |
        # (segs.gstart < BD.x2(boundix) & segs.gend > BD.x2(boundix)) |
        # (segs.gstart >= BD.x1(boundix)& segs.gend <= BD.x2(boundix))
        # ) &
        # ismember(segs.Sample,patient)

        c1 = ((segs['gstart'] <= bd_x1) & (segs['gend'] >= bd_x1))
        c2 = ((segs['gstart'] < bd_x2) & (segs['gend'] > bd_x2))
        c3 = ((segs['gstart'] >= bd_x1) & (segs['gend'] <= bd_x2))
        c4 = (segs['Sample'] == patient)

        segix = (c1 | c2 | c3) & c4
        if sum(segix) == 0:
            continue

        seg1 = segs.loc[segix]

        # sum(seg1.length) < (BD.x2(boundix) - BD.x1(boundix))/arm_length_fraction_threshold

        if (seg1['length'].sum() < ((bd_x2 - bd_x1) / arm_length_fraction_threshold)):
            continue

        region_median = mf.calc_region_median(seg1, bd_x1, bd_x2, 2)

        # region_median < single_del_threshold & region_median >= double_del_threshold

        if (region_median < single_del_threshold) & (region_median >= double_del_threshold):
            arm_del_df.loc[arm + '.DEL', patient] = 1
            #arm_del_df.loc[arm + '.DEL.CCF', patient] = seg1.loc[seg1['Segment_Mean'] == region_median, 'CCF_hat'].max()
        elif region_median < double_del_threshold:
            arm_del_df.loc[arm + '.DEL', patient] = 2
            #arm_del_df.loc[arm + '.DEL.CCF', patient] = seg1.loc[seg1['Segment_Mean'] == region_median, 'CCF_hat'].max()


# Arm amps
arm_amp_df = pd.DataFrame(0, index=BA['arm'] + '.AMP', columns=sorted(segs['Sample'].unique()))
#arm_amp_df_ccf = arm_amp_df.copy(deep=True)
#arm_amp_df_ccf.index = arm_amp_df_ccf.index + '.CCF'
#arm_amp_df = pd.concat([arm_amp_df, arm_amp_df_ccf])

# Use loops for now just to ensure code exact replication.
# This is VERY slow, and should be vectorized.

for patient in arm_amp_df.columns:
    for arm in BA['arm']:

        boundix = (BA['arm'] == arm)
        boundix = boundix[boundix].index

        ba_x1 = BA.loc[boundix, 'x1'].values[0]
        ba_x2 = BA.loc[boundix, 'x2'].values[0]

        # (
        # (segs.gstart <= BD.x1(boundix) & segs.gend >= BD.x1(boundix)) |
        # (segs.gstart < BD.x2(boundix) & segs.gend > BD.x2(boundix)) |
        # (segs.gstart >= BD.x1(boundix)& segs.gend <= BD.x2(boundix))
        # ) &
        # ismember(segs.Sample,patient)

        c1 = ((segs['gstart'] <= ba_x1) & (segs['gend'] >= ba_x1))
        c2 = ((segs['gstart'] < ba_x2) & (segs['gend'] > ba_x2))
        c3 = ((segs['gstart'] >= ba_x1) & (segs['gend'] <= ba_x2))
        c4 = (segs['Sample'] == patient)

        segix = (c1 | c2 | c3) & c4
        if sum(segix) == 0:
            continue

        seg1 = segs.loc[segix]

        # sum(seg1.length) < (BD.x2(boundix) - BD.x1(boundix))/arm_length_fraction_threshold
        if seg1['length'].sum() < ((ba_x2 - ba_x1) / arm_length_fraction_threshold):
            continue

        region_median = mf.calc_region_median(seg1, ba_x1, ba_x2, 2)

        if (region_median > single_amp_threshold) & (region_median <= double_amp_threshold):
            arm_amp_df.loc[arm + '.AMP', patient] = 1
            #arm_amp_df.loc[arm + '.AMP.CCF', patient] = seg1.loc[seg1['Segment_Mean'] == region_median, 'CCF_hat'].max()
        elif region_median > double_amp_threshold:
            arm_amp_df.loc[arm + '.AMP', patient] = 2
            #arm_amp_df.loc[arm + '.AMP.CCF', patient] = seg1.loc[seg1['Segment_Mean'] == region_median, 'CCF_hat'].max()

# Focals
focal_df = pd.DataFrame(0, index=AL['Descriptor'], columns=sorted(segs['Sample'].unique()))
focal_df.insert(0, 'cohort', AL['cohort'].values)
#focal_df_ccf = focal_df.copy(deep=True)
#focal_df_ccf.index = focal_df_ccf.index + '.CCF'
#focal_df = pd.concat([focal_df, focal_df_ccf])

for patient in focal_df.columns[1::]:
    print(patient)
    seg1 = segs.loc[segs['Sample'] == patient]

    if seg1.shape[0] == 0:
        continue

    for peak in AL['Descriptor']:
        if peak[1] in {'p', 'q'}:
            arm = peak[0:2]
        else:
            arm = peak[0:3]

        # (seg1.gstart < AL.gstart(alix) & seg1.gend > AL.gstart(alix)) |
        # (seg1.gstart < AL.gend(alix) & seg1.gend > AL.gend(alix)) |
        # (seg1.gstart > AL.gstart(alix) & seg1.gend < AL.gend(alix));

        gstart_alix = AL.loc[AL['Descriptor'] == peak, 'gstart'].values[0]
        gend_alix = AL.loc[AL['Descriptor'] == peak, 'gend'].values[0]

        c1 = ((seg1['gstart'] < gstart_alix) & (seg1['gend'] > gstart_alix))
        c2 = ((seg1['gstart'] < gend_alix) & (seg1['gend'] > gend_alix))
        c3 = ((seg1['gstart'] > gstart_alix) & (seg1['gend'] < gend_alix))

        segix = c1 | c2 | c3
        segix = segix[segix]

        seg2 = seg1.loc[segix.index]

        if seg2.shape[0] == 0:
            continue

        if ((arm in BA['arm'].values) & ('AMP' in peak)) | ((arm in BD['arm'].values) & ('DEL' in peak)):
            bound1 = arm_level_significance.loc[arm_level_significance['arm'] == arm, 'x1'].values[0]
            bound2 = arm_level_significance.loc[arm_level_significance['arm'] == arm, 'x2'].values[0]
            arm_median = mf.calc_region_median(seg1, bound1, bound2, 2)
            seg2['Segment_Mean'] = seg2['Segment_Mean'] - arm_median

        if 'AMP' in peak:
            seg2['Segment_Mean'] = seg2['Segment_Mean'] * -1

        region_median = mf.calc_region_median(seg2, gstart_alix, gend_alix, 5)

        if 'AMP' in peak:
            region_median = region_median * -1

        # (strfindk(working_matrix.gene(i2), 'AMP') & region_median > double_amp_threshold) |
        # (strfindk(working_matrix.gene(i2), 'DEL') & region_median < double_del_threshold)

        if (('AMP' in peak) and (region_median > double_amp_threshold)) or (('DEL' in peak) and (region_median < double_del_threshold)):
            focal_df.loc[peak, patient] = 2
            """
            # max(seg2.CCF_hat((seg2.Segment_Mean == region_median) | (seg2.Segment_Mean == -region_median)));
            m1 = seg2.loc[seg2['Segment_Mean'] == region_median, 'CCF_hat'].max()
            m2 = seg2.loc[seg2['Segment_Mean'] == -region_median, 'CCF_hat'].max()
            if np.isnan(m1):
                ccf = m2
            elif np.isnan(m2):
                ccf = m1
            else:
                ccf = max(m1, m2)
            focal_df.loc[peak + '.CCF', patient] = ccf
            """
        # (strfindk(working_matrix.gene(i2),'AMP') & (region_median <= double_amp_threshold & region_median > single_amp_threshold)) |
        # (strfindk(working_matrix.gene(i2),'DEL') & (region_median >= double_del_threshold & region_median < single_del_threshold));
        elif (('AMP' in peak) and (region_median <= double_amp_threshold) and (region_median > single_amp_threshold)) or \
             (('DEL' in peak) and (region_median >= double_del_threshold) and (region_median) < single_del_threshold):
            focal_df.loc[peak, patient] = 1
            """
            m1 = seg2.loc[seg2['Segment_Mean'] == region_median, 'CCF_hat'].max()
            m2 = seg2.loc[seg2['Segment_Mean'] == -region_median, 'CCF_hat'].max()
            if np.isnan(m1):
                ccf = m2
            elif np.isnan(m2):
                ccf = m1
            else:
                ccf = max(m1, m2)
            focal_df.loc[peak + '.CCF', patient] = ccf
			"""
focal_df.index = focal_df.index.str.replace(':', '.')
focal_df = focal_df.loc[focal_df['cohort'].str.contains('DLBCL')]

scna_df = pd.concat([arm_del_df, arm_amp_df, focal_df])
scna_df.index = scna_df.index.str.upper()

scna_df.to_csv(OUTPUT_FN, sep='\t')