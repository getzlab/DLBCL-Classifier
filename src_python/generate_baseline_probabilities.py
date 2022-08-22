import pandas as pd
import numpy as np
import scipy.stats as ss
from matplotlib import pyplot as plt

ALPHA = 1
POWER = 2

version = 'Aug_17_2022'

connFile = '../clustering_runs/full_dir_combined/k5/GSM699_cluster_' + version + '.k5.connectivity.matrix.txt'
clusterFile = '../data_tables/clustering_labels/GSM699_cluster_' + version + '.bestclus.remapped.txt'
purityFile = '../data_tables/purities_ploidies/ALLPurities_fixednames.tsv'
coverageFile = '../data_tables/clustering_labels/per_sample_tumor_coverage_formatted.txt'

connMat = pd.read_csv(connFile, delimiter='\t')
clusters = pd.read_csv(clusterFile, delimiter='\t', header=0)
clusters = clusters.sort_values('cluster', axis=0).reset_index(drop=True)
clusters.set_index('SampleName', inplace=True)
purities = pd.read_csv(purityFile, delimiter='\t', header=None, index_col=0)

# clean up the coverages file
coverages = pd.read_csv(coverageFile, delimiter='\t', header=0, index_col = 0)
connMat = connMat.loc[clusters.index.values,clusters.index.values]

c1Subset = clusters.loc[clusters['cluster'] == 1]
c2Subset = clusters.loc[clusters['cluster'] == 2]
c3Subset = clusters.loc[clusters['cluster'] == 3]
c4Subset = clusters.loc[clusters['cluster'] == 4]
c5Subset = clusters.loc[clusters['cluster'] == 5]

probMatrix = pd.DataFrame(0, index=connMat.index.values, columns=['C1', 'C2', 'C3', 'C4', 'C5'])

inconsistent_samples = []
for sample in connMat:
    sortedCol = connMat[sample].sort_values(ascending=False)
    cluster = clusters.loc[sample, 'cluster']

    c1Prob = np.mean(connMat.loc[sample, c1Subset.index.values].values)
    c2Prob = np.mean(connMat.loc[sample, c2Subset.index.values].values)
    c3Prob = np.mean(connMat.loc[sample, c3Subset.index.values].values)
    c4Prob = np.mean(connMat.loc[sample, c4Subset.index.values].values)
    c5Prob = np.mean(connMat.loc[sample, c5Subset.index.values].values)

    probArray = [c1Prob, c2Prob, c3Prob, c4Prob, c5Prob]

    # Fix inconsistencies in connectivity matrix

    if cluster == 1 and np.argmax(probArray) != 0:
        c1Prob = np.max(probArray) + 0.01
        inconsistent_samples.append(sample)

    if cluster == 2 and np.argmax(probArray) != 1:
        c2Prob = np.max(probArray) + 0.01
        inconsistent_samples.append(sample)

    if cluster == 3 and np.argmax(probArray) != 2:
        c3Prob = np.max(probArray) + 0.01
        inconsistent_samples.append(sample)

    if cluster == 4 and np.argmax(probArray) != 3:
        c4Prob = np.max(probArray) + 0.01
        inconsistent_samples.append(sample)

    if cluster == 5 and np.argmax(probArray) != 4:
        c5Prob = np.max(probArray) + 0.01
        inconsistent_samples.append(sample)

    # renormalize
    allSum = c1Prob+c2Prob+c3Prob+c4Prob+c5Prob
    c1Prob = c1Prob/allSum
    c2Prob = c2Prob/allSum
    c3Prob = c3Prob/allSum
    c4Prob = c4Prob/allSum
    c5Prob = c5Prob/allSum

    # Purity penalty term
    purity = purities.loc[sample, 1]
    coverage = coverages.loc[sample, 'coverage_median']
    sensitivity = 1-ss.binom.cdf(3, int(coverage), purity*0.5)
    penaltyTerm = 1-sensitivity
    c1Prob += penaltyTerm
    c2Prob += penaltyTerm
    c3Prob += penaltyTerm
    c4Prob += penaltyTerm
    c5Prob += penaltyTerm

    # Stretch
    c1Prob = c1Prob ** POWER
    c2Prob = c2Prob ** POWER
    c3Prob = c3Prob ** POWER
    c4Prob = c4Prob ** POWER
    c5Prob = c5Prob ** POWER

    # renormalize
    allSum = c1Prob + c2Prob + c3Prob + c4Prob + c5Prob
    c1Prob = c1Prob / allSum
    c2Prob = c2Prob / allSum
    c3Prob = c3Prob / allSum
    c4Prob = c4Prob / allSum
    c5Prob = c5Prob / allSum

    probMatrix.loc[sample, 'C1'] = c1Prob
    probMatrix.loc[sample, 'C2'] = c2Prob
    probMatrix.loc[sample, 'C3'] = c3Prob
    probMatrix.loc[sample, 'C4'] = c4Prob
    probMatrix.loc[sample, 'C5'] = c5Prob

print(inconsistent_samples)
print(len(inconsistent_samples))

probMatrix['cluster'] = probMatrix.idxmax(axis=1).map({'C1': 1, 'C2': 2, 'C3': 3, 'C4': 4, 'C5': 5})
probMatrix['confidence'] = probMatrix.iloc[:, 0:5].max(axis=1)
print(probMatrix['cluster'])

probMatrix = probMatrix.sort_values(by=['cluster', 'confidence'], ascending=[True, False])

probMatrix.to_csv("../data_tables/confidence_tables/baseline_probabilities.connectivity_based.sensitivity_power" +
                  str(POWER) + '.' + version + ".tsv",
                  header=True, index=True, sep="\t")
