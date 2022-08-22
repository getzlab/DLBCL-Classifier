from sklearn.naive_bayes import MultinomialNB
from sklearn.naive_bayes import GaussianNB
import pandas as pd
import generate_test_train_sets as gt
import format_data as fd
import numpy as np
import random


def main(seed, folds, no_sv=False, no_cna=False, ploidy=False, no_silent=False, full_features=False,
         reduced_version=None, use_pca=None, coo=False, subset=False, no_muts=False, qval=None, training_file=None,
         remove_outliers=False):
    # Set seeds
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
        t_targets_arr.append(list(train_targets.idxmax(axis=1)))
        v_inputs_arr.append(validation_data.values)
        v_targets_arr.append(list(validation_targets.idxmax(axis=1)))

    all_preds = []
    all_prob_preds = []
    prob_vectors = []
    samples = []
    models = []
    for k in range(0, folds):
        mnb = MultinomialNB(alpha=1)
        if use_pca:
            mnb = GaussianNB()
        mnbmodel = mnb.fit(t_inputs_arr[k], t_targets_arr[k])
        y_pred = mnbmodel.predict_proba(v_inputs_arr[k])
        for i in range(0, len(y_pred)):
            ele = y_pred[i]
            all_preds.append(ele.argmax())
            all_prob_preds.append(ele.max())
            prob_vectors.append(ele)
            samples.append(validation_sets[k][i])
        models.append(mnbmodel)

    return models, prob_vectors, samples
