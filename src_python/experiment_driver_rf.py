import train_rf as trf
import argparse
import sys
import pickle
import shutil
import os
import pandas as pd

NUM_ITER = 1
folds = 5

parser = argparse.ArgumentParser()
parser.add_argument("--numiter", help='Number of iterations to run', type=int)
parser.add_argument("--folds", help='Number of folds in cross validation', type=int)
parser.add_argument("--nosv", help='Remove all SVs from data table', action='store_true')
parser.add_argument("--nocna", help='Remove all CNAs from data table', action='store_true')
parser.add_argument("--ploidy", help='Include Genome Doubling as a feature', action='store_true')
parser.add_argument("--nosilent", help='Remove all Silent mutations from data table', action='store_true')
parser.add_argument("--binarylabels", help='Use categorical target labels', action='store_true')
parser.add_argument("--fullfeatures", help='Use full (152) feature matrix', action='store_true')
parser.add_argument("--pca", help='Use PCA with n components on feature set before training', type=int)
parser.add_argument("--verbose", help='print more training info', action='store_true')
parser.add_argument("--reduced", help='Use biologically reduced feature set version X', type=str, default=None)
parser.add_argument("--qval", help='Restrict features by q val cutoff', type=float, default=0.10)
parser.add_argument("--coo", help='Use Cell of Origin as a feature.', action='store_true')
parser.add_argument("--subset", help='downsample to samples with COO data', action='store_true')
parser.add_argument("--savemodels", help='Save the models', action='store_true')
parser.add_argument("--nomuts", help='No mutations', action='store_true')
parser.add_argument("--trainingfile", help='Specify a different training set', type=str, default='../data_tables/train_test_sets/TrainingSet_550Subset_May2021.txt')
parser.add_argument("--removeoutliers", help='Remove significantly different (by cohort) markers.', action='store_true')

args = parser.parse_args()

if args.numiter:
    NUM_ITER = args.numiter
if args.folds:
    folds = args.folds

outputFileName = 'RF_evaluation_seeds1_'+str(NUM_ITER) + "_folds"+str(folds)
if args.nosv:
    outputFileName = outputFileName + "_no.sv"
if args.nocna:
    outputFileName = outputFileName + "_no.cna"
if args.ploidy:
    outputFileName = outputFileName + "_ploidy"
if args.nosilent:
    outputFileName = outputFileName + "_no.silent"
if args.binarylabels:
    outputFileName = outputFileName + "_binary.labels"
if args.fullfeatures:
    outputFileName = outputFileName + "_full.features"
if args.pca:
    outputFileName = outputFileName + "_pca" + str(args.pca)
if args.reduced:
    outputFileName = outputFileName + "_reducedV" + str(args.reduced)
if args.qval != 0.10:
    outputFileName = outputFileName + "_qval" + str(args.qval)
if args.coo:
    outputFileName = outputFileName + "_coo"
if args.subset:
    outputFileName = outputFileName + "_subset"
if args.removeoutliers:
    outputFileName = outputFileName + '_remove.outliers'

if args.savemodels:
    dirpath = '../saved_models/' + outputFileName
    if os.path.exists(dirpath) and os.path.isdir(dirpath):
        shutil.rmtree(dirpath)
    os.mkdir(dirpath)


NFEATURES = None
sample_predictions = pd.DataFrame()
for i in range(1, NUM_ITER + 1):
    print("training rf with seed "+str(i))
    forests, validationSets, validationInputs, validationTargets = trf.train_rfs(i, folds, 1,
                                                                                 no_sv=args.nosv,
                                                                                 no_cna=args.nocna,
                                                                                 ploidy=args.ploidy,
                                                                                 no_silent=args.nosilent,
                                                                                 binary_labels=args.binarylabels,
                                                                                 full_features=args.fullfeatures,
                                                                                 use_pca=args.pca,
                                                                                 coo=args.coo,
                                                                                 reduced_version = args.reduced,
                                                                                 subset=args.subset,
                                                                                 no_muts=args.nomuts,
                                                                                 qval=args.qval,
                                                                                 training_file=args.trainingfile,
                                                                                 remove_outliers=args.removeoutliers)

    n_forests = len(forests[0])
    predictions, prob_predictions, classification = trf.predict_rfs(forests, validationInputs)
    if not NFEATURES:
        NFEATURES = len(validationInputs[0][0])

    for pred_idx in range(0, len(predictions)):
        for pred_iter in range(0, len(predictions[pred_idx])):
            pred = predictions[pred_idx][pred_iter]
            sample = validationSets[pred_idx][pred_iter]
            if sample in sample_predictions.columns:
                sample_predictions[sample] += pred
            else:
                sample_predictions[sample] = pred

    if args.savemodels:
        for j in range(0, len(forests)):
            forest = forests[j]
            num = (j+1) + (i-1) * args.folds
            forest_fn = dirpath + '/' + outputFileName + '_' + str(num)
            with open(forest_fn, 'wb+') as f:
                pickle.dump(forest, f, -1)

sample_predictions = sample_predictions.div(NUM_ITER)

outputFileName = outputFileName + '_nfeatures' + str(NFEATURES)
outputFileName = outputFileName + ".tsv"

sample_predictions.to_csv('../evaluation_validation_set/' + outputFileName, sep='\t')
