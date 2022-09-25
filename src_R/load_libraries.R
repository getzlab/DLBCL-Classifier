library(ggplot2)
library(reshape2)
library(tidyr)
library(plyr)
library(dplyr)
library(gtable)
library(stringr)
library(cowplot)
library(VennDiagram)

getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}