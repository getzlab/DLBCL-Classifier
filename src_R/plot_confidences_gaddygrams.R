rm(list = ls())
source('src_R/load_libraries.R')

conf_table = read.csv('evaluation_test_set/NN_reducedV3.4_removeN5_nfeatures21_testsetEval.tsv', 
                      sep='\t', row.names=1)
conf_table = conf_table[order(conf_table$PredictedCluster, conf_table$Confidence), ]
cohorts = read.csv('data_tables/sample_sets/ShippStaudtSets.purity0.2.txt', sep='\t', row.names=1)
cohorts = cohorts[rownames(cohorts) %in% rownames(conf_table),]
cohorts = cohorts[rownames(conf_table), ]

conf_table$cohort = cohorts$cohort
group.colors = c("#FF10F0","#000000")

output_fn = 'plots/test_set/confidences_gaddygram_testset'

conf_table$breaks = seq(1,nrow(conf_table))

p_test <- ggplot(conf_table, aes(x=factor(seq(1,nrow(conf_table))),y=Confidence)) +
  geom_point(aes(shape=Correctness), size=3, alpha=0.5) +
  ggtitle('Cluster-Sorted Confidences (Test Set)') +
  theme_bw() +
  ylim(0,1) +
  xlab('Sample') +
  ylab('Confidence') +
  theme(axis.text.x = element_text(colour="grey20",size=10,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="grey20",size=10,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="grey20",size=20,angle=0,hjust=.5,vjust=0,face="plain"),
        axis.title.y = element_text(colour="grey20",size=20,hjust=.5,vjust=.5,face="plain"),
        plot.title = element_text(hjust=0.50, size=15)) +
  scale_x_discrete(breaks=c(1,25,68,94,118), 
                   labels=c("C1", "C2", "C3", 'C4', 'C5')) +
  expand_limits(x=c(-10, 159)) +
  scale_colour_manual(values = group.colors[1:2]) +
  scale_shape_manual(values = c(17,4)) +
  guides(colour = guide_legend(override.aes = list(shape = 13))) +
  guides(shape = guide_legend(override.aes = list(shape = c(17, 4))))

ggsave(paste(output_fn, '.jpeg', sep=''), p_test, width=8, height=8)
ggsave(paste(output_fn, '.pdf', sep=''), p_test, width=8, height=8)

##############
# Training set
##############

conf_table_train = read.csv('evaluation_validation_set/confidence_adjusted_tables/NN_reducedV3.4_removeN5_nfeatures21_pMax0.93856484.tsv', 
                      sep='\t', row.names=1)
conf_table_train$Confidence = apply(conf_table_train, 1, max)
conf_table_train$PredictedCluster = apply(conf_table_train[,1:5], 1, which.max)
conf_table_train$PredictedCluster = mapvalues(conf_table_train$PredictedCluster,
                                              from=c(1,2,3,4,5),
                                              to=c('C1', 'C2', 'C3', 'C4', 'C5'))
conf_table_train = conf_table_train[order(conf_table_train$PredictedCluster, conf_table_train$Confidence), ]

labels = read.csv('data_tables/confidence_tables/baseline_probabilities.connectivity_based.sensitivity_power2.Sep_23_2022.tsv',
                  sep='\t', row.names=1)
labels = labels[rownames(conf_table_train),]
conf_table_train$TrueCluster = labels$cluster
conf_table_train$TrueCluster = mapvalues(conf_table_train$TrueCluster,
                                         from=c(1,2,3,4,5),
                                         to=c('C1', 'C2', 'C3', 'C4', 'C5'))

conf_table_train$Correctness = (conf_table_train$PredictedCluster == conf_table_train$TrueCluster)
conf_table_train$Correctness = mapvalues(conf_table_train$Correctness, from=c(TRUE, FALSE), to=c('True', 'False'))


cohorts = read.csv('data_tables/sample_sets/ShippStaudtSets.purity0.2.txt', sep='\t', row.names=1)
cohorts = cohorts[rownames(cohorts) %in% rownames(conf_table_train),]
cohorts = cohorts[rownames(conf_table_train), ]

conf_table_train$cohort = cohorts$cohort
conf_table_train$cohort = mapvalues(conf_table_train$cohort, 
                                    from=c('Shipp', 'Staudt', 'TCGA'), 
                                    to=c('Chapuy et al.', 'Schmitz et al.', 'Schmitz et al.'))
group.colors = c("#000000","#E3140F")

conf_table_train$breaks = seq(1,nrow(conf_table_train))

p_train <- ggplot(conf_table_train, aes(x=factor(seq(1,nrow(conf_table_train))),y=Confidence)) +
  geom_point(aes(shape=Correctness), size=3, alpha=0.6) +
  ggtitle('Cluster-Sorted Confidences (Train Set)') +
  theme_bw() +
  ylim(0,1) +
  xlab('Sample') +
  ylab('Confidence') + 
  theme(axis.text.x = element_text(colour="grey20",size=10,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="grey20",size=10,hjust=.5,vjust=.5,face="plain"),
        axis.title.x = element_text(colour="grey20",size=20,angle=0,hjust=.5,vjust=0,face="plain"),
        axis.title.y = element_text(colour="grey20",size=20,hjust=.5,vjust=.5,face="plain"),
        plot.title = element_text(hjust=0.50, size=15)) +
  scale_x_discrete(breaks=c(1,114,240,330,403),
                   labels=c("C1", "C2", "C3", 'C4', 'C5')) +
  expand_limits(x=c(-10, 159)) +
  scale_colour_manual(values = group.colors[1:2]) +
  scale_shape_manual(values = c(17,4)) +
  guides(colour = guide_legend(override.aes = list(shape = 13))) +
  guides(shape = guide_legend(override.aes = list(shape = c(17, 4))))

ggsave(paste('plots/confidences_gaddygram_trainset', '.jpeg', sep=''), p_train, width=8, height=8)
ggsave(paste('plots/paper_figures/confidences_gaddygram_trainset', '.pdf', sep=''), p_train, width=8, height=8)

##############
# All samples
##############

conf_table_all = rbind(conf_table[, c('Confidence', 'PredictedCluster', 'TrueCluster')])
conf_table_all$PredictedCluster = paste('C', conf_table_all$PredictedCluster, sep='')
conf_table_all$TrueCluster = paste('C', conf_table_all$TrueCluster, sep='')
conf_table_all = rbind(conf_table_all, conf_table_train[, c('Confidence', 'PredictedCluster', 'TrueCluster')])
conf_table_all = conf_table_all[order(conf_table_all$PredictedCluster, conf_table_all$Confidence), ]
conf_table_all$Correctness = conf_table_all$PredictedCluster == conf_table_all$TrueCluster

hc_all = conf_table_all[conf_table_all$Confidence > 0.70, ]
acc = sum(hc_all$Correctness) / nrow(hc_all)

p_all <- ggplot(conf_table_all, aes(x=factor(seq(1,nrow(conf_table_all))),y=Confidence, color=Correctness)) +
  geom_point(aes(shape=Correctness, size=Correctness, alpha=Correctness), stroke=1) +
  scale_size_manual(values = c(6, 7)) +
  scale_alpha_manual(values = c(1.0, 0.05)) +
  ggtitle('Cluster-Sorted Confidences (All Samples)') +
  theme_bw() +
  ylim(0,1) +
  xlab('Sample') +
  ylab('Confidence') + 
  theme(axis.text.x = element_text(colour="grey20",size=10,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="grey20",size=10,hjust=.5,vjust=.5,face="plain"),
        axis.title.x = element_text(colour="grey20",size=20,angle=0,hjust=.5,vjust=0,face="plain"),
        axis.title.y = element_text(colour="grey20",size=20,hjust=.5,vjust=.5,face="plain"),
        plot.title = element_text(hjust=0.50, size=15)) +
  scale_x_discrete(breaks=c(1,140,308,421,520),
                   labels=c("C1", "C2", "C3", 'C4', 'C5')) +
  expand_limits(x=c(-10, 159)) +
  scale_colour_manual(values = rev(group.colors[1:2])) +
  scale_shape_manual(values = c(4, 20)) +
  guides(colour = guide_legend(override.aes = list(shape = c(4, 20)))) +
  guides(shape = guide_legend(override.aes = list(shape = c(4, 20)))) +
  geom_hline(yintercept=0.70, color='red') +
  geom_text(aes(670, 0.70, label=0.70, vjust=-0.5)) +
  geom_text(aes(550, 0.25, label='509/524 (97.1%)', vjust=1.3))

ggsave(paste('plots/confidences_gaddygram_allsamples', '.jpeg', sep=''), p_all, width=8, height=8)
ggsave(paste('plots/paper_figures/confidences_gaddygram_allsamples', '.pdf', sep=''), p_all, width=8, height=8)

