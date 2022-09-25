from scipy import stats as ss
from pylab import *
import random


def makeBins(X, Y, windowsize=10, step=10, residpower=2):
    X = np.array(X)
    Y = np.array(Y)
    possibleXs = np.linspace(0, 1, windowsize + 1)
    allWindows = []
    lower = 0
    upper = len(X) % windowsize
    if upper == 0:
        upper = windowsize
    # list segmentation lower:upper is inclusive:noninclusive
    while upper <= len(X):
        currX_window = X[lower:upper]
        meanCurr_X = np.mean(currX_window)
        currY_window = Y[lower:upper]
        meanCurr_Y = np.mean(currY_window)
        currBinSize = len(currX_window)
        minBin = np.abs(np.min(currX_window) - meanCurr_X)
        maxBin = np.abs(np.max(currX_window) - meanCurr_X)

        numcorrect = sum(currY_window)
        numincorrect = len(currY_window) - sum(currY_window)
        alpha = numcorrect + 1
        beta = numincorrect + 1
        probability_aboveX = 1 - ss.beta.cdf(meanCurr_X - 0.1, alpha, beta)
        lowerError = ss.beta.ppf(0.15865, alpha, beta)
        upperError = ss.beta.ppf((1 - 0.15865), alpha, beta)

        currWeight = meanCurr_X * (float(currBinSize) / float(windowsize))

        worstcaseresid = np.max([(np.abs(meanCurr_X - 1) ** residpower), (np.abs(meanCurr_X) ** residpower)])
        bestcaseresid = np.min(np.abs(possibleXs - meanCurr_X))

        allWindows.append((meanCurr_X, meanCurr_Y, probability_aboveX, worstcaseresid, bestcaseresid, currWeight, lowerError, upperError, minBin, maxBin))

        lower = upper
        upper += step
    return allWindows


def xyresiduals_window_weighted_fractionPerWindow(X, Y, name, jiggle=None, seed=1, windowsize=20, step=1, format='pdf',
                                saveFileName=None, computeplots=False, showplots=False, residpower=2, bootstrapwindows=True,
                                xlabel='Confidence', ylabel='Accuracy', y_ticks_top=50, y_ticks_bot=5):
    random.seed(seed)
    X = np.array(X)
    Y = np.array(Y)
    allWindows = makeBins(X, Y, windowsize, step, residpower)

    if bootstrapwindows:
        allWindows = sorted(random.choices(allWindows, k=len(allWindows)))

    windowAverages_X =  np.array([tmp for (tmp, _, _, _, _, _, _, _, _, _) in allWindows])
    windowAverages_Y =  np.array([tmp for (_, tmp, _, _, _, _, _, _, _, _) in allWindows])
    cdfs =              np.array([tmp for (_, _, tmp, _, _, _, _, _, _, _) in allWindows])
    worstcases =        np.array([tmp for (_, _, _, tmp, _, _, _, _, _, _) in allWindows])
    bestcases =         np.array([tmp for (_, _, _, _, tmp, _, _, _, _, _) in allWindows])
    allweights =        np.array([tmp for (_, _, _, _, _, tmp, _, _, _, _) in allWindows])
    lowerErrors =       np.array([tmp for (_, _, _, _, _, _, tmp, _, _, _) in allWindows])
    upperErrors =       np.array([tmp for (_, _, _, _, _, _, _, tmp, _, _) in allWindows])
    lowerBin =          np.array([tmp for (_, _, _, _, _, _, _, _, tmp, _) in allWindows])
    upperBin =          np.array([tmp for (_, _, _, _, _, _, _, _, _, tmp) in allWindows])

    upperErrors = upperErrors - windowAverages_Y
    upperErrors = np.maximum(upperErrors, [0])
    lowerErrors = windowAverages_Y - lowerErrors
    lowerErrors = np.maximum(lowerErrors, [0])

    residuals = (np.abs(windowAverages_X - windowAverages_Y) ** 2)
    window_weightedkappas = (cdfs * (worstcases - residuals) / (worstcases - bestcases)) * allweights
    weighted_meankappa = ((sum((window_weightedkappas)) / sum(allweights)) + 1) / 2

    if computeplots:
        xticks_all = np.linspace(0,1,11)
        fig, axs = plt.subplots(3, gridspec_kw={'height_ratios': [1, 3, 1]})
        for i in range(len(windowAverages_X)):
            x = [windowAverages_X[i]]
            y = [windowAverages_Y[i]]
            xerr = [[lowerBin[i]], [upperBin[i]]]
            yerr = [[lowerErrors[i]], [upperErrors[i]]]
            axs[1].errorbar(x, y, xerr=xerr, yerr=yerr, c='black', ls='none', marker='o', capsize=0, ms=3, alpha=max([allweights[i] - 0.15, 0.20]))
        # axs[1].errorbar(windowAverages_X, windowAverages_Y, xerr=[lowerBin, upperBin], yerr=[lowerErrors, upperErrors],
        #                 c='black', ls='none', marker='o', capsize=0, ms=3)
        axs[1].text(0, 0.80, 'Kappa = ' + str(np.round(weighted_meankappa, 4)), color='blue', size=10)
        axs[1].set_xlim([0, 1])
        axs[1].set_ylim([-0.07, 1.07])
        axs[1].set_ylabel(ylabel)
        axs[1].set_xticks([])
        axs[1].grid(alpha=0.45)
        axs[1].set_axisbelow(True)
        for tick in xticks_all:
            axs[1].axvline(tick, c='gray', alpha=0.45, lw=1, zorder=0)

        axs[0].set_xticks([])
        axs[0].grid(alpha=0.45)
        axs[0].set_axisbelow(True)
        axs[0].set_xlim([0, 1])
        y_h, x_h, _ = axs[0].hist(X[Y.astype(bool)], bins=25, color='green', alpha=0.5, range=(0.0, 1.0))
        upper_y = y_h.max() + round(y_h.max() * 0.30)
        upper_y = int(round(upper_y/5.0)*5.0)
        axs[0].set_ylim([0, upper_y])
        axs[0].set_yticks(np.arange(0, upper_y + 1, y_ticks_top))
        for tick in xticks_all:
            axs[0].axvline(tick, c='gray', alpha=0.45, lw=1, zorder=0)

        #axs[2].set_position([0.125, 0.08, .778, 0.1])
        axs[2].set_xticks([])
        #axs[2].set_ylim(axs[0].get_ylim())
        axs[2].set_xlim([0, 1])
        y_h, x_h, _ = axs[2].hist(X[~Y.astype(bool)], bins=25, color='red', alpha=0.5, range=(0.0, 1.0))
        upper_y = y_h.max() + round(y_h.max() * 0.30)
        upper_y = int(round(upper_y / 5.0) * 5.0)
        axs[2].set_ylim([0, upper_y])
        axs[2].set_yticks(np.arange(0, upper_y + 1, y_ticks_bot))
        axs[2].invert_yaxis()
        axs[2].grid(alpha=0.45)
        axs[2].set_axisbelow(True)
        axs[2].set_xlabel(xlabel)
        axs[2].set_xticks(xticks_all)
        for tick in xticks_all:
            axs[2].axvline(tick, c='gray', alpha=0.45, lw=1, zorder=0)

        x2 = np.linspace(-5, 5, 100)
        axs[1].plot(x2, x2, '-m')
        #plt.title(name)
        if showplots:
            plt.show()
        if saveFileName is not None:
            if saveFileName is not None:
                saveFileName = saveFileName + '.' + format
                plt.savefig(saveFileName, format=format)
        if jiggle is not None:
            Y_jiggle = Y + jiggle
            axs[1].scatter(X[Y.astype(bool)], Y_jiggle[Y.astype(bool)], c='green', s=50, alpha=0.3)
            axs[1].scatter(X[~Y.astype(bool)], Y_jiggle[~Y.astype(bool)], c='red', s=50, alpha=0.3)
        else:
            axs[1].scatter(X[Y.astype(bool)], Y[Y.astype(bool)], c='green', s=50, alpha=0.5)
            axs[1].scatter(X[~Y.astype(bool)], Y[~Y.astype(bool)], c='red', s=50, alpha=0.5)
        if showplots:
            plt.show()
        if saveFileName is not None:
            if saveFileName is not None:
                saveFileName = saveFileName + '_withcorrectness.' + format
                plt.savefig(saveFileName, format=format)
        plt.close()

    return weighted_meankappa, cdfs, residuals, window_weightedkappas
