import pandas as pd
import numpy as np
import correlation_plot as CP
import glob
import re
from operator import itemgetter
import sklearn.metrics
import scipy.stats as ss
import random
from matplotlib import pyplot as plt
from multiprocessing import Pool
import os

# Reduce the # of pool workers that get allocated...
# This parameter is necessary to avoid MAC OS crashes

# default values:
# BOOTSTRAP_SAMPLES = True
# DEFAULT_SEED = 123
# NUM_CORRELATION_ITERS = 1000
# WINDOW_SIZE = 10
# STEP_SIZE = 10


DEFAULT_SEED = 123
NUM_CORRELATION_ITERS = 1000
WINDOW_SIZE = 10
STEP_SIZE = 10

random.seed(DEFAULT_SEED)
np.random.seed(DEFAULT_SEED)

# list files in ensemble validation evaluation
files = glob.glob('../evaluation_validation_set/*')
files = list(reversed(sorted(files)))


# Make regex for files in evaluation_validation_set
nameReg = re.compile('.+(NN|RF|MNB)_evaluation_seeds\d_\d+_folds\d_(.+).tsv')
nameReg2 = re.compile('.+(NN|RF|MNB)_evaluation_noretrain_(.+).tsv')


def evaluate_model(file):
    random.seed(DEFAULT_SEED)
    np.random.seed(DEFAULT_SEED)
    # Make regex for files in EnsembleValidationEvaluation
    targetfile = "../data_tables/confidence_tables/baseline_probabilities.connectivity_based.sensitivity_power2.Aug_17_2022.tsv"
    labels = pd.read_csv(targetfile, sep='\t', index_col=0)
    labels = pd.DataFrame(labels['cluster'] - 1)

    name_reg = re.compile('.+(NN|RF|MNB)_evaluation_seeds\d_\d+_folds\d_(.+).tsv')
    name_reg2 = re.compile('.+(NN|RF|MNB)_evaluation_noretrain_(.+).tsv')
    tuples = []
    name = name_reg.findall(file)
    name2 = name_reg2.findall(file)
    if not name:
        if not name2:
            return ''
        else:
            name = name2

    name = '%s_%s' % (name[0][0], name[0][1])
    name = name.replace('reducedV4', 'qval0.10')

    print('Evaluating', name)
    predictions_df = pd.read_csv(file, sep='\t', index_col=0)
    predictions_df = predictions_df.T
    # do confidence transform
    p_max = np.max(predictions_df.max())
    for idx, row in predictions_df.iterrows():
        new_arr = np.float64(np.array(row))
        p_win = np.float64(np.max(row))
        new_conf = np.float64(p_win / p_max)
        shrinkfactor = np.float64((p_max - p_win) / (p_max * np.float64(np.sum(np.delete(new_arr, new_arr.argmax())))))
        for idx2 in range(len(new_arr)):
            if idx2 == new_arr.argmax():
                new_arr[idx2] = new_conf
            else:
                new_arr[idx2] = new_arr[idx2] * shrinkfactor
        predictions_df.loc[idx] = new_arr
    predictions_df.to_csv('../evaluation_validation_set/confidence_adjusted_tables/' + name + '_pMax' + str(p_max) + '.tsv', sep='\t')

    for idx, row in predictions_df.iterrows():
        tuples.append((idx, row.max(), row.idxmax(), labels.loc[idx]['cluster']))
    tuples = sorted(tuples, key=itemgetter(1))
    samples = np.array([x for (x, _, _, _) in tuples])
    confidences = np.array([x for (_, x, _, _) in tuples])
    predclus = np.array([x for (_, _, x, _) in tuples])
    trueclus = np.array([x for (_, _, _, x) in tuples])

    top70th_actuals = trueclus[int(len(trueclus) * .30)::]
    top70th_preds = predclus[int(len(predclus) * .30)::]
    top50th_actuals = trueclus[int(len(trueclus) * .50)::]
    top50th_preds = predclus[int(len(predclus) * .50)::]
    top15th_actuals = trueclus[int(len(trueclus) * .85)::]
    top15th_preds = predclus[int(len(predclus) * .85)::]
    top5th_actuals = trueclus[int(len(trueclus) * .95)::]
    top5th_preds = predclus[int(len(predclus) * .95)::]

    above90conf = confidences > 0.90
    above90conf_preds = predclus[above90conf]
    above90conf_actuals = trueclus[above90conf]
    b_80_90conf = np.equal(confidences <= 0.90, confidences > 0.80)
    b_80_90conf_preds = predclus[b_80_90conf]
    b_80_90conf_actuals = trueclus[b_80_90conf]
    b_70_80conf = np.equal(confidences <= 0.80, confidences > 0.70)
    b_70_80conf_preds = predclus[b_70_80conf]
    b_70_80conf_actuals = trueclus[b_70_80conf]
    below70conf = confidences <= 0.70
    below70conf_preds = predclus[below70conf]
    below70conf_actuals = trueclus[below70conf]

    confusion_matrix_top = sklearn.metrics.confusion_matrix(top70th_actuals, top70th_preds)
    accuracy_top70 = np.trace(confusion_matrix_top)
    accuracy_top70 = accuracy_top70 / np.sum(confusion_matrix_top)
    confusion_matrix_top50 = sklearn.metrics.confusion_matrix(top50th_actuals, top50th_preds)
    accuracy_top50 = np.trace(confusion_matrix_top50)
    accuracy_top50 = accuracy_top50 / np.sum(confusion_matrix_top50)
    confusion_matrix_top15 = sklearn.metrics.confusion_matrix(top15th_actuals, top15th_preds)
    accuracy_top15 = np.trace(confusion_matrix_top15)
    accuracy_top15 = accuracy_top15 / np.sum(confusion_matrix_top15)
    confusion_matrix_top5 = sklearn.metrics.confusion_matrix(top5th_actuals, top5th_preds)
    accuracy_top5 = np.trace(confusion_matrix_top5)
    accuracy_top5 = accuracy_top5 / np.sum(confusion_matrix_top5)
    all_accuracy = sum(np.equal(trueclus, predclus)) / len(predclus)
    confusion_matrix_all = sklearn.metrics.confusion_matrix(trueclus, predclus)

    confusion_matrix_above90 = sklearn.metrics.confusion_matrix(above90conf_actuals, above90conf_preds)
    accuracy_above90 = np.trace(confusion_matrix_above90) / np.sum(confusion_matrix_above90)

    confusion_matrix_b_80_90 = sklearn.metrics.confusion_matrix(b_80_90conf_actuals, b_80_90conf_preds)
    accuracy_b_80_90 = np.trace(confusion_matrix_b_80_90) / np.sum(confusion_matrix_b_80_90)

    confusion_matrix_b_70_80 = sklearn.metrics.confusion_matrix(b_70_80conf_actuals, b_70_80conf_preds)
    accuracy_b_70_80 = np.trace(confusion_matrix_b_70_80) / np.sum(confusion_matrix_b_70_80)

    confusion_matrix_below70 = sklearn.metrics.confusion_matrix(below70conf_actuals, below70conf_preds)
    accuracy_below70 = np.trace(confusion_matrix_below70) / np.sum(confusion_matrix_below70)

    with open('../data_tables/confusion_matrices/' + name + '.txt', 'w+') as confF:
        confF.write('Format: TRUE labels in ROWS, PREDICTED labels in COLUMNS\n')
        confF.write('Confusion Matrix All\n')
        writestr = np.array2string(confusion_matrix_all, separator='\t')
        writestr = writestr.replace('[', ' ')
        writestr = writestr.replace(']', ' ')
        confF.write(writestr)
        confF.write('\nOverall Accuracy: ' + str(all_accuracy))
        confF.write('\n')
        confF.write('Confusion Matrix Top 70%\n')
        writestr = np.array2string(confusion_matrix_top, separator='\t')
        writestr = writestr.replace('[', ' ')
        writestr = writestr.replace(']', ' ')
        confF.write(writestr)
        confF.write('\nAccuracy Top 70%: ' + str(accuracy_top70))
        confF.write('\n')
        confF.write('Confusion Matrix Above 0.90\n')
        writestr = np.array2string(confusion_matrix_above90, separator='\t')
        writestr = writestr.replace('[', ' ')
        writestr = writestr.replace(']', ' ')
        confF.write(writestr)
        confF.write('\nAccuracy Above 0.90: ' + str(accuracy_above90))
        confF.write('\n')
        confF.write('Confusion Matrix 0.80_0.90\n')
        writestr = np.array2string(confusion_matrix_b_80_90, separator='\t')
        writestr = writestr.replace('[', ' ')
        writestr = writestr.replace(']', ' ')
        confF.write(writestr)
        confF.write('\nAccuracy 0.80 < Conf <= 0.90: ' + str(accuracy_b_80_90))
        confF.write('\n')
        confF.write('Confusion Matrix 0.70_0.80\n')
        writestr = np.array2string(confusion_matrix_b_70_80, separator='\t')
        writestr = writestr.replace('[', ' ')
        writestr = writestr.replace(']', ' ')
        confF.write(writestr)
        confF.write('\nAccuracy 0.70 < Conf <= 0.80: ' + str(accuracy_b_70_80))
        confF.write('\n')
        confF.write('Confusion Matrix <= 0.70\n')
        writestr = np.array2string(confusion_matrix_below70, separator='\t')
        writestr = writestr.replace('[', ' ')
        writestr = writestr.replace(']', ' ')
        confF.write(writestr)
        confF.write('\nAccuracy Conf <= 0.70: ' + str(accuracy_below70))

    correctness = np.equal(trueclus, predclus)
    zipped = list(zip(correctness, confidences, samples))
    lrx = np.array([tup[1] for tup in zipped])
    lry = np.array([tup[0] for tup in zipped], dtype=int)

    # Compute CI on top overall accuracy
    alpha = sum(correctness) + 1
    beta = len(correctness) - (alpha - 1) + 1
    lower_acc = ss.beta.ppf(0.15865, alpha, beta)
    upper_acc = ss.beta.ppf((1 - 0.15865), alpha, beta)

    mean_top_acc = np.mean([accuracy_top50, accuracy_top70, accuracy_top15, accuracy_top5])

    _, _, _, _ = CP.xyresiduals_window_weighted_fractionPerWindow(lrx, lry, name, seed=1, bootstrapwindows=False,
                                                                  windowsize=WINDOW_SIZE, step=STEP_SIZE, computeplots=True,
                                                                  saveFileName='../plots/xyresiduals/' + name,
                                                                  format='pdf')

    meankappa_val, cdfs_val, residuals_val, weightedkappas_val = \
        CP.xyresiduals_window_weighted_fractionPerWindow(lrx, lry, name, seed=1, bootstrapwindows=False,
                                                         windowsize=WINDOW_SIZE, step=STEP_SIZE, computeplots=True,
                                                         saveFileName='../plots/xyresiduals/' + name,
                                                         format='png')
    ########### y = x residual evaluation ############
    kappas = []
    mean_resids = []
    for i in range(NUM_CORRELATION_ITERS):
        meankappa_boot, cdfs_boot, residuals_boot, weightedkappas_boot = \
            CP.xyresiduals_window_weighted_fractionPerWindow(lrx, lry, name, seed=i,
                                                             windowsize=WINDOW_SIZE, step=STEP_SIZE, bootstrapwindows=True)
        kappas.append(meankappa_boot)
        mean_resids.append(residuals_boot)

    meanresidual = np.mean(mean_resids)
    corr = np.mean(kappas)

    plt.clf()
    plt.hist(kappas, bins=100)
    plt.axvline(x=meankappa_val, color='r')
    plt.title(name + "_kappas_bootstrapped")
    plt.xlabel('Bootstrapped Kappas')
    plt.savefig('../plots/xyresiduals/bootstrap_distributions/' + name + "_kappas_hist.png", format='png')
    plt.savefig('../plots/xyresiduals/bootstrap_distributions/' + name + "_kappas_hist.pdf", format='pdf')
    plt.close()

    BETA = 2
    lower_corr = corr - np.std(kappas)
    upper_corr = corr + np.std(kappas)
    perf_numer = all_accuracy * corr
    perf_denom = all_accuracy + (corr * (BETA ** 2))

    # Compute F1/F_beta score
    # error on accuracy (binomial error) is asymmetric so compute upper/lower errors separately

    perf_numer_err_lower = np.sqrt(((corr - lower_corr) / corr) ** 2 +
                                ((all_accuracy - lower_acc) / all_accuracy) ** 2)
    perf_numer_err_upper = np.sqrt(((upper_corr - corr) / corr) ** 2 +
                                ((upper_acc - all_accuracy) / all_accuracy) ** 2)

    perf_denom_err_lower = np.sqrt((corr - lower_corr) ** 2 + (all_accuracy - lower_acc) ** 2)
    perf_denom_err_upper = np.sqrt((upper_corr - corr) ** 2 + (upper_acc - all_accuracy) ** 2)

    # division/multiplication error propagation uses % error instead of absolute error, so we need to
    # multiply the computed error percentages by the actual value itself
    perf_numer_err_lower = perf_numer * perf_numer_err_lower
    perf_numer_err_upper = perf_numer * perf_numer_err_upper
    perf_denom_err_lower = perf_denom_err_lower
    perf_denom_err_upper = perf_denom_err_upper

    perf_err_upper_percent = np.sqrt((perf_numer_err_upper / perf_numer) ** 2 + (perf_denom_err_upper / perf_denom) ** 2)
    perf_err_lower_percent = np.sqrt((perf_numer_err_lower / perf_numer) ** 2 +
                                  (perf_denom_err_lower / perf_denom) ** 2)

    performance = (1 + BETA ** 2) * perf_numer / perf_denom
    perf_err_upper = perf_err_upper_percent * performance
    perf_err_lower = perf_err_lower_percent * performance
    upper_perf = performance + perf_err_upper
    lower_perf = performance - perf_err_lower

    return ('\n' + name + '\t' +
            str(all_accuracy) + '\t' + str(accuracy_top70) + '\t' + str(accuracy_top50) + '\t' +
            str(accuracy_top15) + '\t' + str(accuracy_top5) + '\t' +
            str(meanresidual) + '\t' + str(corr) + '\t' + str(performance) + '\t' +
            str(lower_acc) + '\t' + str(upper_acc) + '\t' +
            str(lower_corr) + '\t' + str(upper_corr) + '\t' +
            str(lower_perf) + '\t' + str(upper_perf))


with Pool(len(files)) as p:
    outputs = p.map(evaluate_model, files)


with open('../evaluation_validation_set/allmodelsevaluated.tsv', 'w+') as f:
    f.write('experiment\t' +
            'accuracyAll\taccuracyTop70percent\taccuracyTop50percent\taccuracyTop15percent\taccuracyTop5percent\t' +
            'mean_residual\tKappa\tperformance\t' +
            'lowerAcc\tupperAcc\t' +
            'lowerKappa\tupperKappa\t' +
            'lowerPerformance\tupperPerformance')
    for output in outputs:
        if output:
            f.write(output)
