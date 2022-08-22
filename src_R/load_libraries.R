library(ggplot2)
library(reshape2)
library(tidyr)
library(plyr)
library(dplyr)
library(gtable)
library(stringr)
library(cowplot)

getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}




# library(grid)
# library(randomForest)
# library(caret)
# library(ggplot2)
# library(reshape2)
# library(tidyr)
# library(plyr)
# library(dplyr)
# library(stats)
# library(Hmisc)
# library(Rtsne)
# library(readxl)
# library(gridExtra)
# library(gtable)
# library(plotrix)
# library(tsne)
# library(circlize)
# library(ggforce)
# library(stringr)
# 
# # if (!requireNamespace("BiocManager", quietly = TRUE))
# # install.packages("BiocManager")
# 
# # BiocManager::install("DNAcopy")
# # library(DNAcopy)
# 
# library(cluster)
# 
# # if (!requireNamespace("BiocManager", quietly = TRUE))
# #install.packages("BiocManager")
# 
# #BiocManager::install("pd.genomewidesnp.6")
# library(pd.genomewidesnp.6)
# 
# library(magrittr)
# library(jsonlite)
# library(Biostrings)
# library(betareg)
# library(cowplot)
# library(ggrepel)
# library(ggpubr)
# 
# getmode <- function(v) {
#   uniqv <- unique(v)
#   uniqv[which.max(tabulate(match(v, uniqv)))]
# }
# 
# mapclusters <- function(labels, matrix){
#   colMax <- function(data) sapply(data, which.max)
#   hmatMaxes = colMax(matrix)
#   clusters = as.vector(as.numeric(labels$cluster[!duplicated(labels$cluster)]))
#   mapvec = integer(length(clusters))
#   for(cluster in clusters){
#     clusterSamples = rownames(labels)[labels$cluster == cluster]
#     mappedcluster = getmode(as.vector(hmatMaxes[clusterSamples]))
#     mapvec[cluster] = mappedcluster
#   }
#   return(mapvec)
# }