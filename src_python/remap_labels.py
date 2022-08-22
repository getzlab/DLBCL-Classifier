import sys
import pandas as pd
import re
import scipy.stats as ss

filename = sys.argv[1]
outputfilename = re.split('\.\w+$', filename)
outputfilename = outputfilename[0]
outputfilename = outputfilename + '.remapped.txt'

originalLabels = pd.read_csv('../data_tables/clustering_labels/NatMed.DLBCL.bestclus.txt', delimiter='\t')
newLabels = pd.read_csv(filename, delimiter='\t', skiprows=1)
newLabels['unmapped_cluster'] = newLabels['cluster']

clusters = list(set(newLabels['cluster']))

if len(clusters) != 5:
    print('WARNING: Number of clusters is not 5!! Will continue to map clusters.')

clusterMap = {}
for clus in clusters:
    clusSamples = newLabels.loc[newLabels['cluster'] == clus, 'SampleName'].values
    originalVals = originalLabels.loc[originalLabels.index.isin(clusSamples), 'cluster'].values
    mappedCluster = ss.mode(originalVals)[0][0]
    clusterMap[clus] = mappedCluster

newLabels['cluster'] = newLabels['cluster'].map(clusterMap).values
newLabels = newLabels.sort_values(by='cluster')
newLabels.to_csv(outputfilename, sep='\t', header=True, index=False)
print("Wrote file to: "+outputfilename)