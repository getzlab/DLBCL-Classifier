import train_nn as tnn
import argparse
import sys
import pandas as pd
import numpy as np
from operator import itemgetter
import shutil
import os
import torch

NUM_ITER = 1
folds = 5

parser = argparse.ArgumentParser()
parser.add_argument("--numiter", help='Number of iterations to run', type=int)
parser.add_argument("--folds", help='Number of folds in cross validation', type=int)
parser.add_argument("--trimfeatures", help='Trim the least significant n features', type=int)
parser.add_argument("--nosv", help='Remove all SVs from data table', action='store_true')
parser.add_argument("--nocna", help='Remove all CNAs from data table', action='store_true')
parser.add_argument("--ploidy", help='Include Genome Doubling as a feature', action='store_true')
parser.add_argument("--nosilent", help='Remove all Silent mutations from data table', action='store_true')
parser.add_argument("--binarylabels", help='Use categorical target labels', action='store_true')
parser.add_argument("--fullfeatures", help='Use full feature matrix', action='store_true')
parser.add_argument("--pca", help='Use PCA with n components on feature set before training', type=int)
parser.add_argument("--verbose", help='print more training info', action='store_true')
parser.add_argument("--earlystopping", help='each Epoch, check validation error.', action='store_true', default=True)
parser.add_argument("--reduced", help='Use biologically reduced feature set version X', type=str, default=None)
parser.add_argument("--qval", help='Restrict features by q val cutoff', type=float, default=0.10)
parser.add_argument("--coo", help='Use Cell of Origin as a feature.', action='store_true')
parser.add_argument("--subset", help='downsample to samples with COO data', action='store_true')
parser.add_argument("--savemodels", help='Save the models', action='store_true')
parser.add_argument("--nomuts", help='No mutations', action='store_true')
parser.add_argument("--noarms", help='No arm level events', action='store_true')
parser.add_argument("--nofocals", help='No focal level events', action='store_true')
parser.add_argument("--bpcutoff", help='Cut features by cumulative sum <= bpcutoff', type=int, default=None)
parser.add_argument("--keepSVs", help='Keep SVs instead of dropping (for bpcutoff only)', action='store_true')
parser.add_argument("--traininghistory", help='Store the training/validation loss history', action='store_true')
parser.add_argument("--trainingfile", help='Specify a different training set', type=str, default='../data_tables/train_test_sets/TrainingSet_550Subset_May2021.txt')
parser.add_argument("--validationsets", help='Generate validation sets. Does not train model.', action='store_true')
parser.add_argument("--removeoutliers", help='Remove significantly different (by cohort) markers.', action='store_true')
parser.add_argument("--weightl265p", help='Weight MYD88 L265P to 3 instead of 2 for nonsilent.', action='store_true')
parser.add_argument("--remove_largest_n", help='Remove the largest (by footprint) n features.', type=int, default=None)
parser.add_argument("--nosvbcl6", help='Remove the largest (by footprint) n features.', action='store_true')


args = parser.parse_args()

if args.numiter:
    NUM_ITER = args.numiter
if args.folds:
    folds = args.folds

outputFileName = 'NN_evaluation_seeds1_'+str(NUM_ITER) + "_folds"+str(folds)
if args.nosv:
    outputFileName = outputFileName + "_no.sv"
if args.nocna:
    outputFileName = outputFileName + "_no.cna"
if args.nomuts:
    outputFileName = outputFileName + "_no.muts"
if args.noarms:
    outputFileName = outputFileName + "_no.arms"
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
if args.bpcutoff:
    outputFileName = outputFileName + "_bp.cutoff" + str(args.bpcutoff)
if args.keepSVs:
    outputFileName = outputFileName + "_keep.svs"
if args.qval != 0.10:
    outputFileName = outputFileName + "_qval" + str(args.qval)
if args.coo:
    outputFileName = outputFileName + "_coo"
if args.subset:
    outputFileName = outputFileName + "_subset"
if args.nofocals:
    outputFileName = outputFileName + "_no.focals"
if args.trimfeatures:
    outputFileName = outputFileName + "_trim.features" + str(args.trimfeatures)
if args.trainingfile != '../data_tables/train_test_sets/TrainingSet_550Subset_May2021.txt':
    outputFileName = outputFileName + "_pseudo.test.set"
if args.removeoutliers:
    outputFileName = outputFileName + '_remove.outliers'
if args.weightl265p:
    outputFileName = outputFileName + '_weightl265p'
if args.remove_largest_n:
    outputFileName = outputFileName + '_removeN' + str(args.remove_largest_n)
if args.nosvbcl6:
    outputFileName = outputFileName + '_nosvBCL6'


performances, pearsonrs, accuraciesTop = [], [], []
n_features = None
predictionDF = pd.DataFrame()

if args.savemodels:
    savepath = '../saved_models/' + outputFileName
    if os.path.exists(savepath) and os.path.isdir(savepath):
        shutil.rmtree(savepath)
    os.mkdir(savepath)

if args.traininghistory:
    dirpath = '../model_training_history/' + outputFileName
    if os.path.exists(dirpath) and os.path.isdir(dirpath):
        shutil.rmtree(dirpath)
    os.mkdir(dirpath)

if args.validationsets:
    valpath = '../all_validation_sets/' + outputFileName
    if os.path.exists(valpath) and os.path.isdir(valpath):
        shutil.rmtree(valpath)
    os.mkdir(valpath)

samplePredictions = None

for i in range(1, NUM_ITER + 1):
    print("Running train_nn.py with seed "+str(i))
    nets, actualConfidences, predConfidences, allpredictions = \
        tnn.main(i, args.folds, no_sv=args.nosv, no_cna=args.nocna, ploidy=args.ploidy, no_silent=args.nosilent, binary_labels=args.binarylabels,
                 full_features=args.fullfeatures, reduced_version=args.reduced, verbose=args.verbose, early_stopping=args.earlystopping,
                 use_pca=args.pca, trim_features=args.trimfeatures, coo=args.coo, subset=args.subset,
                 no_muts=args.nomuts, generate_validation_sets=args.validationsets, output_filename=outputFileName,
                 training_history=args.traininghistory, qval=args.qval, no_arms=args.noarms, no_focals=args.nofocals, bp_cutoff=args.bpcutoff,
                 keep_svs=args.keepSVs, training_file=args.trainingfile, remove_outliers=args.removeoutliers, weightl265p=args.weightl265p,
                 remove_largest_n=args.remove_largest_n, nosvbcl6=args.nosvbcl6)
    if args.validationsets:
        continue
    if not n_features:
        n_features = nets[0].inputLayer.in_features

    if samplePredictions is None:
        samplePredictions = pd.DataFrame()
        for ele in allpredictions:
            samplePredictions[ele[4]] = ele[5][0]
    else:
        for ele in allpredictions:
            samplePredictions[ele[4]] += ele[5][0]

    if args.savemodels:
        for j in range(0, len(nets)):
            net = nets[j]
            num = (j+1) + (i-1) * args.folds
            netfn = savepath + '/' + outputFileName + '_' + str(num)
            torch.save(net.state_dict(), netfn)

    allpredictions = sorted(allpredictions, key=itemgetter(0))
    c1confidences = [x for (x, _, y, _, _, _) in allpredictions if y == 0]
    c2confidences = [x for (x, _, y, _, _, _) in allpredictions if y == 1]
    c3confidences = [x for (x, _, y, _, _, _) in allpredictions if y == 2]
    c4confidences = [x for (x, _, y, _, _, _) in allpredictions if y == 3]
    c5confidences = [x for (x, _, y, _, _, _) in allpredictions if y == 4]
    fullconfidences = [c1confidences, c2confidences, c3confidences, c4confidences, c5confidences]
    clusters = [[1 for _ in c1confidences],
                [2 for _ in c2confidences],
                [3 for _ in c3confidences],
                [4 for _ in c4confidences],
                [5 for _ in c5confidences]]

    for ele in allpredictions:
        sample = ele[4]
        outputvec = ele[5][0]
        outputvec = np.append(outputvec, ele[2])
        outputvec = np.append(outputvec, ele[3])
        if sample not in predictionDF.index:
            entry = pd.DataFrame([outputvec], columns=['C1', 'C2', 'C3', 'C4', 'C5',
                                                       'ActualCluster', 'PredictedCluster'], index=[sample])
            predictionDF = predictionDF.append(entry)
        else:
            predictionDF.loc[sample] += outputvec


if args.validationsets:
    exit(0)

samplePredictions = samplePredictions.div(NUM_ITER)

outputFileName = outputFileName + '_nfeatures' + str(n_features)
outputFileName = outputFileName + ".tsv"

samplePredictions.to_csv('../evaluation_validation_set/'+outputFileName, sep='\t')
