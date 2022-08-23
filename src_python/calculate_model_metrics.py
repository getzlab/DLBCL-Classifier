import numpy as np
import scipy.stats as ss
import correlation_plot as CP
import random

def modelAccuracy(predcluster, truecluster):
    correctness = np.equal(predcluster, truecluster)
    alpha = sum(correctness) + 1
    beta = len(correctness) - (alpha - 1) + 1
    accuracy = sum(correctness) / len(correctness)

    # edge cases for acc is 0 or 1

    if accuracy == 0:
        lowerAcc = 0
    else:
        lowerAcc = ss.beta.ppf(0.15865, alpha, beta)

    if accuracy == 1:
        upperAcc = 1
    else:
        upperAcc = ss.beta.ppf((1 - 0.15865), alpha, beta)

    return accuracy, lowerAcc, upperAcc


def modelKappa(confidences, correctness, num_iter=1000, name='NoName', window_size=10, step_size=10, bootstrapwindows=True, seed=123,
               computeplots=False, showplots=False, savefilename=None, format=None):
    np.random.seed(seed)
    random.seed(seed)
    kappas = []
    # fractions_bootstrapped = []
    for i in range(num_iter):
        meankappa_boot, cdfs_boot, residuals_boot, weightedkappas_boot = \
            CP.xyresiduals_window_weighted_fractionPerWindow(confidences, correctness, name, seed=i,
                                                             windowsize=window_size, step=step_size, bootstrapwindows=bootstrapwindows,
                                                             computeplots=computeplots, showplots=showplots, saveFileName=savefilename)
        kappas.append(meankappa_boot)

    Kappa = np.mean(kappas)
    lowerKappa = Kappa - np.std(kappas)
    upperKappa = Kappa + np.std(kappas)
    return Kappa, lowerKappa, upperKappa


def modelPerformance(accuracy, kappa, lowerAccuracy, upperAccuracy, lowerKappa, upperKappa):
    BETA = 2
    perfNumer = accuracy * kappa
    perfDenom = accuracy + (kappa * (BETA ** 2))

    # Compute F1/F_beta score
    # error on accuracy (binomial error) is asymmetric so compute upper/lower errors separately
    # BETA = 1

    perfNumerErrLower = np.sqrt(((kappa - lowerKappa) / kappa) ** 2 +
                                ((accuracy - lowerAccuracy) / accuracy) ** 2)
    perfNumerErrUpper = np.sqrt(((upperKappa - kappa) / kappa) ** 2 +
                                ((upperAccuracy - accuracy) / accuracy) ** 2)

    perfDenomErrLower = np.sqrt((kappa - lowerKappa) ** 2 + (accuracy - lowerAccuracy) ** 2)
    perfDenomErrUpper = np.sqrt((upperKappa - kappa) ** 2 + (upperAccuracy - accuracy) ** 2)

    # division/multiplication error propagation uses % error instead of absolute error, so we need to
    # multiply the computed error percentages by the actual value itself
    perfNumerErrLower = perfNumer * perfNumerErrLower
    perfNumerErrUpper = perfNumer * perfNumerErrUpper
    perfDenomErrLower = perfDenomErrLower
    perfDenomErrUpper = perfDenomErrUpper

    perfErrUpperPercent = np.sqrt((perfNumerErrUpper / perfNumer) ** 2 + (perfDenomErrUpper / perfDenom) ** 2)
    perfErrLowerPercent = np.sqrt((perfNumerErrLower / perfNumer) ** 2 +
                                  (perfDenomErrLower / perfDenom) ** 2)

    performance = (1 + BETA ** 2) * perfNumer / perfDenom
    perfErrUpper = perfErrUpperPercent * performance
    perfErrLower = perfErrLowerPercent * performance
    upperPerf = performance + perfErrUpper
    lowerPerf = performance - perfErrLower

    return performance, lowerPerf, upperPerf






