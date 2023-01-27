rm(list = ls())
source("src_R/load_libraries.R")

labels = read.csv('data_tables/confidence_tables/baseline_probabilities.connectivity_based.sensitivity_power2.Sep_23_2022.tsv',
                  sep='\t', row.names=1)

full_gsm = read.csv('data_tables/gsm/DLBCL.699.fullGSM.Sep_23_2022.tsv', sep='\t', row.names=1)

ploidy = full_gsm['PLOIDY', rownames(labels)]
ploidy['Sample', ] = colnames(ploidy)
ploidy = data.frame(t(ploidy))
ploidy$PLOIDY = as.double(ploidy$PLOIDY)
ploidy$cluster = paste('C', labels$cluster, sep='')

p <- ggplot(ploidy, aes(x=cluster, y=PLOIDY)) + 
  geom_jitter(height=0, width=0.15) +
  theme() + theme_bw() + ylab('Ploidy')


plotheight = 4
plotwidth = 8

ggsave('plots/paper_figures/cluster_ploidy.pdf', p, height = plotheight, width = plotwidth)
ggsave('plots/cluster_ploidy.png', p, height = plotheight, width = plotwidth)
