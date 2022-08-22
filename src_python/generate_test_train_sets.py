import pandas as pd
import numpy as np
import random


def generate_kfold_frames(targetfile, k, current_samples, seed):
    np.random.seed(seed)
    random.seed(seed)
    labels = pd.read_csv(targetfile, delimiter='\t', index_col=0)
    labels = labels[labels.index.isin(current_samples)]
    labels = labels.sample(frac=1)

    clusters = list(set(labels['cluster']))
    nclusters = len(clusters)

    samples = [[]]*nclusters
    for c in clusters:
        cluster_subset = labels[labels['cluster'] == c].index.values
        idx = c-1
        samples[idx] = cluster_subset

    validation_sets = [[] for _ in range(k)]
    training_sets = [[] for _ in range(k)]
    for clus in samples:
        fold = 0
        for s in clus:
            foldidx = fold % k
            validation_sets[foldidx].append(s)
            fold += 1

    all_samples = set(labels.index.values)
    for i in range(len(validation_sets)):
        curr_train_set = list(all_samples-set(validation_sets[i]))
        training_sets[i] = curr_train_set

    return training_sets, validation_sets


def generate_train_validation_frames(data_frame, target_frame, training_sets, validation_sets, k):
    training_samples, validation_samples = training_sets[k], validation_sets[k]

    train_data = data_frame[data_frame.index.isin(training_samples)]
    train_targets = target_frame[target_frame.index.isin(training_samples)]

    train_targets = train_targets.reindex(train_data.index.values)

    validation_data = data_frame[data_frame.index.isin(validation_samples)]
    validation_targets = target_frame[target_frame.index.isin(validation_samples)]

    validation_data = validation_data.reindex(validation_samples)
    validation_targets = validation_targets.reindex(validation_data.index.values)

    return train_data, train_targets, validation_data, validation_targets


def generate_test_frames(data_frame, target_frame, test_set):
    test_data = data_frame.loc[test_set]
    test_targets = target_frame.loc[test_set]

    return test_data, test_targets
