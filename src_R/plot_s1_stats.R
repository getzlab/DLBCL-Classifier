rm(list = ls())
source('src_R/load_libraries.R')

s1_table = read.csv('data_tables/tableS1_classifier_merged.tsv', sep='\t', row.names = 1)
full_set = read.csv('data_tables/confidence_tables/baseline_probabilities.connectivity_based.sensitivity_power2.Sep_23_2022.tsv', 
                    sep='\t', row.names=1)
gsm = read.csv('data_tables/gsm/DLBCL.699.fullGSM.Sep_23_2022.tsv',
               sep='\t', row.names=1)
drivers = read.csv('data_tables/qval_dfs/fisher_exact_5x2.Sep_23_2022.combined.tsv', sep='\t', row.names=1)

raw_counts = read.csv('data_tables/raw_muts_cnas_counts.tsv', sep='\t', row.names=1)
raw_sv = read.csv('data_tables/raw_sv_counts.tsv', sep='\t', row.names=1)

pvals_clinical = read.csv('data_tables/tableS2A_clinical.tsv', sep='\t', row.names=1)
pvals_genomic = read.csv('data_tables/tableS2B_genomic.tsv', sep='\t', row.names=1)

driver_set = rownames(drivers[drivers$q < 0.10, ])
bl = c("OR51B6", "OR10V1")
driver_set = driver_set[!driver_set %in% bl]


gsm_drivers = mutate_all(gsm[driver_set, ], function(x) as.numeric(as.character(x)))
amps = driver_set[grepl('\\.AMP', driver_set)]
dels = driver_set[grepl('\\.DEL', driver_set)]
cnas = c(amps, dels)
svs = driver_set[grepl('SV\\.', driver_set)]
muts = driver_set[(!driver_set %in% cnas) & (!driver_set %in% svs)]

# sets = read.csv('data_tables/sample_sets/ShippStaudtSets.purity0.2.txt', sep='\t', row.names=1)
# 
# sets_AS = sets[rownames(sets) %in% rownames(full_set),]
# sets_AS$cohort2 = ifelse(sets_AS$cohort == 'Shipp', 'Chapuy et al.', 'Schmitz et al.')

s1_table_AS = s1_table[rownames(s1_table) %in% rownames(full_set), ]

purity_meds = aggregate(s1_table_AS$Purity..ABSOLUTE., list(s1_table_AS$Cohort), FUN=median)
ploidy_meds = aggregate(s1_table_AS$Ploidy..ABSOLUTE., list(s1_table_AS$Cohort), FUN=median)
cov_meds = aggregate(s1_table_AS$Mean.Target.Coverage, list(s1_table_AS$Cohort), FUN=median)

colnames(purity_meds) = c('Cohort', 'Median_Purity')
colnames(ploidy_meds) = c('Cohort', 'Median_Ploidy')
colnames(cov_meds) = c('Cohort', 'Median_Mean_Target_Coverage')

purity_meds$SD = aggregate(s1_table_AS$Purity..ABSOLUTE., list(s1_table_AS$Cohort), FUN=sd)$x
ploidy_meds$SD = aggregate(s1_table_AS$Ploidy..ABSOLUTE., list(s1_table_AS$Cohort), FUN=sd)$x
cov_meds$SD = aggregate(s1_table_AS$Mean.Target.Coverage, list(s1_table_AS$Cohort), FUN=sd)$x

plot_pur = ggplot(purity_meds, aes(x=Cohort, y=Median_Purity)) + 
    geom_point() +
    geom_errorbar(aes(ymin=Median_Purity-SD, ymax=Median_Purity+SD), width=.2) +
    scale_y_continuous(limits=c(0.2, 1)) +
    theme_bw()

plot_plo = ggplot(ploidy_meds, aes(x=Cohort, y=Median_Ploidy)) + 
  geom_point() +
  geom_errorbar(aes(ymin=Median_Ploidy-SD, ymax=Median_Ploidy+SD), width=.2) +
  scale_y_continuous(limits=c(1, 3)) +
  theme_bw()


plot_mtc = ggplot(cov_meds, aes(x=Cohort, y=Median_Mean_Target_Coverage)) + 
  geom_point() +
  geom_errorbar(aes(ymin=Median_Mean_Target_Coverage-SD, ymax=Median_Mean_Target_Coverage+SD), width=.2) +
  scale_y_continuous() +
  theme_bw()

comb = plot_grid(plot_pur, plot_plo, plot_mtc, rel_widths = c(1,1,1), ncol = 3)

ggsave('plots/paper_figures/pur_plo_cov.pdf', comb, height = 3, width = 8)
ggsave('plots/pur_plo_cov.png', comb, height = 3, width = 8)

hist_pur = ggplot(s1_table_AS, aes(x=Purity..ABSOLUTE., fill=Cohort)) +
  geom_histogram(bins=30, center=0, alpha=0.7, position = "identity", show.legend = FALSE) +
  scale_x_continuous(limits=c(0, 1)) +
  theme(axis.title.x = element_text()) +
  xlab('Purity') +
  ylab('Number of Tumors') +
  scale_fill_manual(values=c('#3677e0', '#d62728')) +
  theme_bw() +
  geom_text(x=0.15, y=32,
            label=paste('p =', pvals_genomic['Purity (Median)', 'p.value..Schmitz.vs.Chapuy.']),
            color='blue', size=3)


hist_plo = ggplot(s1_table_AS, aes(x=Ploidy..ABSOLUTE., fill=Cohort)) +
  geom_histogram(bins=30, center=0, alpha=0.7, position = "identity", show.legend = FALSE) +
  scale_x_continuous() +
  theme(axis.title.x = element_text(),
        legend.position="none") +
  xlab('Ploidy') +
  scale_fill_manual(values=c('#3677e0', '#d62728')) +
  ylab('') +
  theme_bw() +
  geom_text(x=3.6, y=163.5,
            label=paste('p =', pvals_genomic['Ploidy (Median)', 'p.value..Schmitz.vs.Chapuy.']),
            color='blue', size=3)

hist_mtc = ggplot(s1_table_AS, aes(x=Mean.Target.Coverage, fill=Cohort)) +
  geom_histogram(bins=60, center=0, alpha=0.7, position = "identity") +
  scale_x_continuous() +
  theme(axis.title.x = element_text()) +
  xlab('Mean Target Coverage') +
  scale_fill_manual(values=c('#3677e0', '#d62728')) +
  ylab('') +
  theme_bw() +
  geom_text(x=200, y=31,
            label=paste('p =', pvals_genomic['Mean Target Coverage (Median)', 'p.value..Schmitz.vs.Chapuy.']),
            color='blue', size=3)


comb = plot_grid(hist_pur, hist_plo, hist_mtc, rel_widths = c(1,1,1.7), ncol = 3)

ggsave('plots/paper_figures/hist_pur_plo_cov.pdf', comb, height = 3, width = 8)
ggsave('plots/hist_pur_plo_cov.png', comb, height = 3, width = 8)

###############
# driver sums #
###############

tumor_muts_sums = data.frame(colSums(gsm_drivers[muts,] != 0))
tumor_muts_sums = tumor_muts_sums[rownames(s1_table_AS),,drop=F]
tumor_muts_sums$Cohort = s1_table_AS$Cohort
colnames(tumor_muts_sums) = c('driver_count', 'Cohort')

hist_muts = ggplot(tumor_muts_sums, aes(x=driver_count, fill=Cohort)) +
  geom_histogram(bins=25, alpha=0.7, center=0, show.legend = FALSE, position = "identity") +
  scale_x_continuous() +
  theme(axis.title.x = element_text()) +
  scale_fill_manual(values=c('#3677e0', '#d62728')) +
  xlab('Mutations') +
  ylab('') +
  theme_bw() +
  geom_text(x=33, y=60,
            label=paste('p =', pvals_genomic['Driver Count (Mutations, GSM)', 'p.value..Schmitz.vs.Chapuy.']),
            color='blue', size=3)


tumor_dels_sums = data.frame(colSums(gsm_drivers[dels,] != 0))
tumor_dels_sums = tumor_dels_sums[rownames(s1_table_AS),,drop=F]
tumor_dels_sums$Cohort = s1_table_AS$Cohort
colnames(tumor_dels_sums) = c('driver_count', 'Cohort')

hist_dels = ggplot(tumor_dels_sums, aes(x=driver_count, fill=Cohort)) +
  geom_histogram(bins=20, alpha=0.7, center=0, show.legend = FALSE, position = "identity") +
  scale_x_continuous() +
  theme(axis.title.x = element_text()) +
  scale_fill_manual(values=c('#3677e0', '#d62728')) +
  xlab('Copy Number Deletions') +
  ylab('') +
  theme_bw() +
  geom_text(x=17, y=70,
            label=paste('p =', pvals_genomic['Driver Count (Deletions, GSM)', 'p.value..Schmitz.vs.Chapuy.']),
            color='blue', size=3)
  
tumor_amps_sums = data.frame(colSums(gsm_drivers[amps, ] != 0))
tumor_amps_sums = tumor_amps_sums[rownames(s1_table_AS),,drop=F]
tumor_amps_sums$Cohort = s1_table_AS$Cohort
colnames(tumor_amps_sums) = c('driver_count', 'Cohort')

hist_amps = ggplot(tumor_amps_sums, aes(x=driver_count, fill=Cohort)) +
  geom_histogram(bins=18, alpha=0.7, center=0, position = "identity") +
  scale_x_continuous() +
  theme(axis.title.x = element_text()) +
  scale_fill_manual(values=c('#3677e0', '#d62728')) +
  xlab('Copy Number Amplifications') +
  ylab('') +
  theme_bw() +
  geom_text(x=14.5, y=47,
            label=paste('p =', pvals_genomic['Driver Count (Amplifications, GSM)', 'p.value..Schmitz.vs.Chapuy.']),
            color='blue', size=3)

tumor_svs_sums = data.frame(colSums(gsm_drivers[svs, ] != 0))
tumor_svs_sums = tumor_svs_sums[rownames(s1_table_AS),,drop=F]
tumor_svs_sums$Cohort = s1_table_AS$Cohort
colnames(tumor_svs_sums) = c('driver_count', 'Cohort')

hist_svs = ggplot(tumor_svs_sums, aes(x=driver_count, fill=Cohort)) +
  geom_histogram(bins=18, alpha=0.7, center=0, show.legend = FALSE, position = "identity") +
  scale_x_continuous() +
  theme(axis.title.x = element_text()) +
  scale_fill_manual(values=c('#3677e0', '#d62728')) +
  xlab('Structual Variants') +
  ylab('Number of Tumors') +
  theme_bw() +
  geom_text(x=2.5, y=250,
            label=paste('p =', pvals_genomic['Driver Count (SVs, GSM)', 'p.value..Schmitz.vs.Chapuy.']),
            color='blue', size=3)


comb = plot_grid(hist_svs, hist_muts, hist_dels, hist_amps, rel_widths = c(1,1,1,1.7), ncol = 4)

ggsave('plots/paper_figures/hist_driver_counts.pdf', comb, height = 3, width = 12)
ggsave('plots/hist_driver_counts.png', comb, height = 3, width = 12)

##############
# raw counts #
##############
raw_counts = raw_counts[rownames(s1_table_AS), ]
raw_counts$Cohort = s1_table_AS$Cohort

muts_c = ggplot(raw_counts, aes(x=muts, fill=Cohort)) +
  geom_histogram(bins=45, center=0, alpha=0.7, show.legend=FALSE, position = "identity") +
  theme(axis.title.x = element_text()) +
  xlab('Mutation Count') +
  ylab('Number of Tumors') +
  scale_fill_manual(values=c('#3677e0', '#d62728')) +
  theme_bw() +
  scale_x_continuous(trans='log10', limits=c(1, 11000)) +
  geom_text(x=0.6, y=63,
            label=paste('p =', pvals_genomic['Number of Mutations (Median)', 'p.value..Schmitz.vs.Chapuy.']),
            color='blue', size=3)

dels_c = ggplot(raw_counts, aes(x=del_frac, fill=Cohort)) +
  geom_histogram(bins=45, center=0, alpha=0.7, show.legend=FALSE, position = "identity") +
  scale_x_continuous() +
  theme(axis.title.x = element_text()) +
  xlab('Fraction Deleted') +
  ylab('') +
  scale_fill_manual(values=c('#3677e0', '#d62728')) +
  theme_bw() +
  geom_text(x=0.5, y=69,
            label=paste('p =', pvals_genomic['Fraction Genome Deleted (Median)', 'p.value..Schmitz.vs.Chapuy.']),
            color='blue', size=3)

amps_c = ggplot(raw_counts, aes(x=amp_frac, fill=Cohort)) +
  geom_histogram(bins=45, center=0, alpha=0.7, position = "identity") +
  scale_x_continuous() +
  theme(axis.title.x = element_text()) +
  scale_fill_manual(values=c('#3677e0', '#d62728')) +
  xlab('Fraction Amplified') +
  ylab('') +
  theme_bw() +
  geom_text(x=0.6, y=35,
            label=paste('p =', pvals_genomic['Fraction Genome Amplified (Median)', 'p.value..Schmitz.vs.Chapuy.']),
            color='blue', size=3)

raw_sv = raw_sv[rownames(s1_table_AS),,drop=FALSE]
raw_sv$Cohort = s1_table_AS$Cohort

svs_c = ggplot(raw_sv, aes(x=individual, fill=Cohort)) +
  geom_histogram(bins=20, alpha=0.7, center=0, show.legend=FALSE) +
  scale_x_continuous() +
  theme(axis.title.x = element_text()) +
  scale_fill_manual(values=c('#3677e0', '#d62728')) +
  xlab('Structural Variations*') +
  ylab('') +
  theme_bw()

comb = plot_grid(svs_c, muts_c, dels_c, amps_c, rel_widths = c(1,1,1,1.6), ncol = 4)

ggsave('plots/paper_figures/hist_raw_counts.pdf', comb, height = 3, width = 12)
ggsave('plots/hist_raw_counts.png', comb, height = 3, width = 12)

