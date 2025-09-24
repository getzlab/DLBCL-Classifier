MAF=$1
SEG=$2
BP=$3


#
# Example file paths assume that this is run from the git directory DLBCL-Classifier/src_python/gsm_generation_scripts
# tested 24Sep2025 using python 3.12 with import pandas, numpy, and argparse 
#
SETID="DLBclass_699"
SETSAMPLE="../../data_tables/train_test_sets/DLBclass_699_Sep2025.txt"
BP="../../data_tables/additional_gsm_inputs/DLBCL_Shipp_Staudt.SVs.14-Dec-2021.txt"
MAF="../../data_tables/maf_files/DLBCL_combined_699.hg38B.noPDE4DIP.noHISTartifacts.trim.maf"
SEG="../../data_tables/additional_gsm_inputs/DLBCL_699_segs.2021-12-15.ccf.seg"
CNVBLACKLIST="../../data_tables/additional_gsm_inputs/CNV.hg19.bypos.111213.CR1_event_added.bed"
GISTICARMS="../../data_tables/additional_gsm_inputs/DLBCL_cHL_PMBL_broad_significance.txt"
GISTICFOCALS="../../data_tables/additional_gsm_inputs/all_significant_focal_peaks_with_cohort.txt"
FEATUREORDER="../../data_tables/gsm/feature_order.19Aug2024X.txt"
OUTDIR="/tmp/"
TODAY=$(date +%d%b%Y)
echo "python3 sv2gsm.py -i $SETID -s $SETSAMPLE -v $BP -o $OUTDIR "
python3 sv2gsm.py --i $SETID -s $SETSAMPLE -v $BP -o $OUTDIR

echo "python3 maf2gsm.py --i $SETID -s $SETSAMPLE -m $MAF -o $OUTDIR"
python3 maf2gsm.py -i $SETID -s $SETSAMPLE -m $MAF -o $OUTDIR

echo "python3 seg2gsm.py -i $SETID -s $SETSAMPLE -v $SEG -x $CNVBLACKLIST -a $GISTICARMS -f $GISTICFOCALS -o $OUTDIR"
python3 seg2gsm.py -i $SETID -s $SETSAMPLE -v $SEG -x $CNVBLACKLIST -a $GISTICARMS -f $GISTICFOCALS -o $OUTDIR


#SVGSM="/tmp/freeze_Aug14_all_248.21Aug2024.SV.GSM.tsv"
#MAFGSM="/tmp/freeze_Aug14_all_248.21Aug2024.MAF.GSM.tsv"
#CNVGSM="/tmp/freeze_Aug14_all_248.21Aug2024.CNV.GSM.tsv"


SVGSM="${OUTDIR}${SETID}.${TODAY}.SV.GSM.tsv" # Concatenates without a space
echo "SV GSM: ${SVGSM}"
MAFGSM="${OUTDIR}${SETID}.${TODAY}.MAF.GSM.tsv" # Concatenates without a space
echo "MAF GSM: ${MAFGSM}"
CNVGSM="${OUTDIR}${SETID}.${TODAY}.CNV.GSM.tsv" # Concatenates without a space
echo "CNVGSM: ${CNVGSM}"



python3 combine2gsm.py -i $SETID -v $SEG -v $SVGSM -m $MAFGSM -c $CNVGSM -f $FEATUREORDER -o $OUTDIR
