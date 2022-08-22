import generate_test_train_sets as gt
import format_data as fd
import numpy as np
import random
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.ensemble import RandomForestRegressor


def train_rfs(seed, folds, nforests,
              no_sv=False, no_cna=False, ploidy=False, no_silent=False, binary_labels=False, full_features=False,
              use_pca=None, coo=False, reduced_version=None, subset=False, no_muts=False, qval=0.10, training_file=None,
              remove_outliers=False):
    np.random.seed(seed)
    random.seed(seed)

    datafile = "../data_tables/gsm/DLBCL_Staudt_Shipp_CL.for_classifier_training.classifier_subset.fix_sv.fix_ploidy.17-Aug-2022.txt"
    targetfile = "../data_tables/confidence_tables/baseline_probabilities.connectivity_based.sensitivity_power2.Aug_17_2022.tsv"

    train_samples = list(pd.read_csv(training_file, sep='\t', header=None, index_col=0).index)

    data, targets = fd.format_inputs(datafile, targetfile, train_samples,
                                     no_sv=no_sv,
                                     no_cna=no_cna,
                                     ploidy=ploidy,
                                     no_silent=no_silent,
                                     reduced_version=reduced_version,
                                     use_pca=use_pca,
                                     coo=coo,
                                     no_muts=no_muts,
                                     qval=qval,
                                     remove_outliers=remove_outliers)

    training_sets, validation_sets = gt.generate_kfold_frames(targetfile, folds, list(data.index), seed)

    t_inputs_arr = []
    t_targets_arr = []
    v_inputs_arr = []
    v_targets_arr = []

    for k in range(0, folds):
        train_data, train_targets, validation_data, validation_targets = \
            gt.generate_train_validation_frames(data, targets, training_sets, validation_sets, k)
        t_inputs_arr.append(train_data.values)
        t_targets_arr.append(train_targets.values)
        v_inputs_arr.append(validation_data.values)
        v_targets_arr.append(validation_targets.values)
        print("Training size: " + str(train_data.shape[0]))

    print("Training on " + str(data.shape[1]) + ' features')

    all_forests = []
    for k in range(0, folds):
        print("Fold: " + str(k+1))
        forests = []
        for f in range(0, nforests):
            curr_labels = t_targets_arr[k]
            labels = []
            if binary_labels:
                for val in curr_labels:
                    labels.append(np.argmax(val))
                rf = RandomForestClassifier(n_estimators=500)
                rf.fit(t_inputs_arr[k], labels)
                forests.append(rf)
            else:
                rf = RandomForestRegressor(n_estimators=500)
                rf.fit(t_inputs_arr[k], curr_labels)
                forests.append(rf)
        all_forests.append(forests)
    return all_forests, validation_sets, v_inputs_arr, v_targets_arr


def predict_rfs(forests, validation_inputs):
    classification = (type(forests[0][0]).__name__ == "RandomForestClassifier")
    predictions = []
    if classification:
        prob_predictions = []
    else:
        prob_predictions = None
    for fold in range(len(validation_inputs)):
        for f in forests[fold]:
            prediction = f.predict(validation_inputs[fold])
            predictions.append(prediction)
            if classification:
                probpred = f.predict_proba(validation_inputs[fold])
                prob_predictions.append(probpred)
    return predictions, prob_predictions, classification
