rm(list = ls())

source('src_R/load_libraries.R')

fullDF = read.csv('data_tables/gsm/DLBCL.699.fullGSM.Sep_23_2022.tsv', 
                  sep='\t', row.names = 1)
qvalDF = read.csv('data_tables/qval_dfs/fisher_exact_5x2.Sep_23_2022.combined.tsv', sep='\t', row.names=1)

qvalDF = qvalDF[!rownames(qvalDF) %in% c("OR10V1" ,"OR51B6"),]
qvalDF = qvalDF[qvalDF$q <= 0.10,]

cohorts = read.csv('data_tables/sample_sets/sample_inclusion_table.tsv', sep='\t', row.names=1)

train_set = read.csv('data_tables/train_test_sets/TrainingSet_550Subset_May2021.txt', sep='\t', header=FALSE)
train_set = train_set$V1

fullDF = data.frame(t(fullDF))
fullDF = fullDF[, !grepl('.CCF', colnames(fullDF))]
fullDF = fullDF[, ]
fullDF = fullDF[, !colnames(fullDF) %in% c('PLOIDY', 'EBV', 'COO')]
cohorts = cohorts[cohorts$included_in_clustering == 'True',]

staudt_coo_downsampled = read.csv('data_tables/phenotypes/staudt_coo_downsample.tsv', sep='\t', row.names=1)

for(col in colnames(fullDF)){
  fullDF[, col] = as.numeric(fullDF[, col])
}

NCI_DF = fullDF[rownames(cohorts)[cohorts$cohort == 'Staudt'],]
NCI_DF = NCI_DF[rownames(NCI_DF) %in% rownames(staudt_coo_downsampled),]
TCGA_DF = fullDF[rownames(fullDF) %in% rownames(cohorts)[cohorts$cohort == 'TCGA'],]

Staudt_DF = rbind(NCI_DF, TCGA_DF)
Our_DF = fullDF[rownames(fullDF) %in% rownames(cohorts)[cohorts$cohort == 'Shipp'], ]

plotDF = data.frame(cbind(colSums(Staudt_DF != 0) / nrow(Staudt_DF), colSums(Our_DF != 0) / nrow(Our_DF)))
colnames(plotDF) = c('Schmitz_F', 'Chapuy_F')
plotDF$diff = plotDF$Schmitz_F - plotDF$Chapuy_F
plotDF$absDiff = abs(plotDF$Schmitz_F - plotDF$Chapuy_F)

plotDF$counts_wt_schmitz = colSums(Staudt_DF == 0)
plotDF$counts_var_schmitz = colSums(Staudt_DF != 0)
plotDF$counts_wt_chapuy = colSums(Our_DF == 0)
plotDF$counts_var_chapuy = colSums(Our_DF != 0)

# Plot DF is just the main dataframe, it contains frequencies for ours/staudt,
# and the counts of variant/wild type events for ours/staudt.
plotDF$pval_fisher = -1

plotDF$AlterationType = 'SNP'
plotDF[grepl('SV', rownames(plotDF)), 'AlterationType'] = 'SV'
plotDF[grepl('\\.AMP', rownames(plotDF)), 'AlterationType'] = 'SCNA - AMP'
plotDF[grepl('\\.DEL', rownames(plotDF)), 'AlterationType'] = 'SCNA - DEL'

for(i in rownames(plotDF)){
  
  fisher_mat = matrix(0, 2, 2)
  fisher_mat[1,1] = plotDF[i, 'counts_wt_schmitz']
  fisher_mat[1,2] = plotDF[i, 'counts_var_schmitz']
  fisher_mat[2,1] = plotDF[i, 'counts_wt_chapuy']
  fisher_mat[2,2] = plotDF[i, 'counts_var_chapuy']
  
  pval_fisher = fisher.test(fisher_mat, alternative = 'two.sided')[1]
  
  plotDF[i, 'pval_fisher'] = pval_fisher
  if (i == 'X21Q22.3.AMP'){
    print(fisher_mat)
  }
}

plotDF$X = rownames(plotDF)

natmedgenes = read.csv('data_tables/gsm/DLBCL_significant_event_matrix_NatMed.txt',
                       sep='\t', row.names = 1)
rownames(natmedgenes)[97] = "X18Q21.32.AMP"
#rownames(natmedgenes)[131] = "X19Q13.32.1.DEL"


rownames(natmedgenes) = toupper(make.names(rownames(natmedgenes)))

plotDF_natmed = plotDF[rownames(plotDF) %in% rownames(natmedgenes),]
plotDF_natmed = plotDF_natmed[!rownames(plotDF_natmed) %in% c("HIST1H2BK", "HIST2H2BE"), ]
plotDF_natmed$qval_natmed_fisher = p.adjust(plotDF_natmed$pval_fisher, method='BH')
plotDF_natmed$significant_fisher = plotDF_natmed$qval_natmed_fisher <= 0.10
plotDF_natmed$logq_fisher = -log10(plotDF_natmed$qval_natmed_fisher)

plotDF = plotDF[rownames(plotDF) %in% rownames(qvalDF),]
plotDF$qval_fisher = p.adjust(plotDF$pval_fisher, method='BH')
plotDF$logq_fisher = -log10(plotDF$qval_fisher)
plotDF$significant_fisher = plotDF$qval_fisher <= 0.10

group.colors <- c("#000000","#E3140F")
alteration.colors = c("#eb004e", "#1cd9d6", "#000000", "#008c59")

write.table(data.frame("Gene"=rownames(plotDF),plotDF),
            "driver_frequencies/frequency_df_staudt_ours.tsv", row.names=FALSE, sep='\t', quote=FALSE)

p <- ggplot(plotDF, aes(x=Chapuy_F, y=Schmitz_F, color=AlterationType, shape=significant_fisher)) + 
  geom_point() +
  xlim(0, 0.4) +
  ylim(0, 0.4) +
  theme_bw() +
  ggtitle(paste('Cohort Frequencies - All Drivers', '(N=', nrow(plotDF) ,')', sep='')) +
  ylab('Schmitz et al. Driver Frequencies') +
  xlab('Chapuy et al. Driver Frequencies') +
  theme(axis.text.x = element_text(colour="grey20",size=10,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="grey20",size=10,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="grey20",size=15,angle=0,hjust=.5,vjust=0,face="plain"),
        axis.title.y = element_text(colour="grey20",size=15,hjust=.5,vjust=.5,face="plain")) +
  geom_text(aes(label=ifelse(qval_fisher <= 0.10, as.character(X),'')), vjust=-0.4, size=2) +
  geom_abline() +
  #geom_abline(slope=1, intercept=0.10, linetype='dashed', color='black', show.legend=TRUE) +
  #geom_abline(slope=1, intercept=-0.10, linetype='dashed', show.legend = TRUE) +
  geom_text(x=0.2, y=0.4, 
            label=paste("Fraction significant (qval <= 0.10) = ", 
                        round(sum(plotDF$qval_fisher <= 0.10) / nrow(plotDF), 3), sep=''),
            color='blue') +
  geom_text(x=0.2, y=0.4-0.03, 
            label=paste('pearson correlation:', round(as.numeric(cor.test(plotDF$Chapuy_F, plotDF$Schmitz_F)[4]), 3), '\n',
                        'p =', format(as.numeric(cor.test(plotDF$Chapuy_F, plotDF$Schmitz_F)[3]), scientific=TRUE)),
            color='blue') +
  scale_colour_manual(values = alteration.colors) +
  scale_shape_manual(values = c(19, 4))

p1 <- ggplot(plotDF_natmed, aes(x=Chapuy_F, y=Schmitz_F, color=AlterationType, shape=significant_fisher)) + 
  geom_point() +
  xlim(0, 0.4) +
  ylim(0, 0.4) +
  theme_bw() +
  ggtitle(paste('Cohort Frequencies - Nat Med Drivers', '(N=', nrow(plotDF_natmed) ,')', sep='')) +
  ylab('Schmitz et al. Driver Frequencies') +
  xlab('Chapuy et al. Driver Frequencies') +
  theme(axis.text.x = element_text(colour="grey20",size=10,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="grey20",size=10,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="grey20",size=15,angle=0,hjust=.5,vjust=0,face="plain"),
        axis.title.y = element_text(colour="grey20",size=15,hjust=.5,vjust=.5,face="plain")) +
  geom_text(aes(label=ifelse(qval_natmed_fisher <= 0.10, as.character(X),'')), vjust=-0.4, size=2) +
  geom_abline() +
  geom_abline(slope=1, intercept=0.10, linetype='dashed') +
  geom_abline(slope=1, intercept=-0.10, linetype='dashed') +
  geom_text(x=0.2, y=0.4, 
            label=paste("Fraction significant (qval <= 0.10) = ", 
                        round(sum(plotDF_natmed$qval_natmed_fisher <= 0.10) / nrow(plotDF_natmed), 3), sep=''),
            color='blue') +
  geom_text(x=0.2, y=0.4-0.03, 
            label=paste('pearson correlation:', round(as.numeric(cor.test(plotDF_natmed$Chapuy_F, plotDF_natmed$Schmitz_F)[4]), 3), '\n',
                        'p =', format(as.numeric(cor.test(plotDF_natmed$Chapuy_F, plotDF_natmed$Schmitz_F)[3]), scientific=TRUE)),
            color='blue') +
  scale_colour_manual(values = alteration.colors) +
  scale_shape_manual(values = c(19, 4))

p2 <- ggplot(plotDF, aes(x=diff, y=logq_fisher, color=AlterationType, shape=significant_fisher)) + 
  geom_point() +
  xlim(-0.2, 0.2) +
  scale_y_continuous(limits=c(0, 8)) +
  theme_bw() +
  ggtitle('Difference vs q value (All Drivers)') +
  xlab('Difference (Staudt - Chapuy)') +
  ylab('-log10(q) : FDR adjusted') +
  geom_text(aes(label=ifelse(qval_fisher <= 0.10, as.character(X),'')), vjust=-0.7, size=2) +
  theme(axis.text.x = element_text(colour="grey20",size=10,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="grey20",size=10,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="grey20",size=15,angle=0,hjust=.5,vjust=0,face="plain"),
        axis.title.y = element_text(colour="grey20",size=15,hjust=.5,vjust=.5,face="plain"),
        plot.title = element_text(hjust=0.50, size=15)) +
  scale_colour_manual(values = alteration.colors) +
  scale_shape_manual(values = c(19, 4))

p3 <- ggplot(plotDF_natmed, aes(x=diff, y=logq_fisher, color=AlterationType, shape=significant_fisher)) + 
  geom_point() +
  xlim(-0.2, 0.2) +
  scale_y_continuous(limits=c(0, 8)) +
  theme_bw() +
  ggtitle('Difference vs q value (Nat. Med. Drivers)') +
  xlab('Difference (Staudt - Chapuy)') +
  ylab('-log10(q) : FDR adjusted') +
  geom_text(aes(label=ifelse(qval_natmed_fisher <= 0.10, as.character(X),'')), vjust=-0.7, size=2) +
  theme(axis.text.x = element_text(colour="grey20",size=10,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="grey20",size=10,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="grey20",size=15,angle=0,hjust=.5,vjust=0,face="plain"),
        axis.title.y = element_text(colour="grey20",size=15,hjust=.5,vjust=.5,face="plain"),
        plot.title = element_text(hjust=0.50, size=15)) +
  scale_colour_manual(values = alteration.colors) +
  scale_shape_manual(values = c(19, 4))

### Train vs Test frequencies

fullDF = fullDF[, colnames(fullDF) %in% rownames(qvalDF)]
fullDF_train_cl = data.frame(t(fullDF[rownames(fullDF) %in% train_set, ]))
fullDF_test_cl = data.frame(t(fullDF[!rownames(fullDF) %in% train_set, ]))


plotDF = data.frame(cbind(rowSums(fullDF_train_cl != 0) / ncol(fullDF_train_cl), 
                          rowSums(fullDF_test_cl != 0) / ncol(fullDF_test_cl)))

colnames(plotDF) = c('Train_F', 'Test_F')

plotDF$counts_wt_train = rowSums(fullDF_train_cl == 0)
plotDF$counts_wt_test = rowSums(fullDF_test_cl == 0)

plotDF$counts_var_train = rowSums(fullDF_train_cl != 0)
plotDF$counts_var_test = rowSums(fullDF_test_cl != 0)

plotDF$pval_fisher = -1

for(i in rownames(plotDF)){
  
  fisher_mat = matrix(0, 2, 2)
  fisher_mat[1,1] = plotDF[i, 'counts_wt_train']
  fisher_mat[1,2] = plotDF[i, 'counts_var_train']
  fisher_mat[2,1] = plotDF[i, 'counts_wt_test']
  fisher_mat[2,2] = plotDF[i, 'counts_var_test']
  
  pval_fisher = fisher.test(fisher_mat, alternative = 'two.sided')[1]
  
  plotDF[i, 'pval_fisher'] = pval_fisher
}

plotDF$qval_fisher = p.adjust(plotDF$pval_fisher, method='BH')
plotDF$logq_fisher = -log10(plotDF$qval_fisher)
plotDF$significant_fisher = plotDF$qval_fisher <= 0.10
plotDF$X = rownames(plotDF)
plotDF$absDiff = abs(plotDF$Train_F - plotDF$Test_F)

plotDF$AlterationType = 'SNP'
plotDF[grepl('SV', rownames(plotDF)), 'AlterationType'] = 'SV'
plotDF[grepl('\\.AMP', rownames(plotDF)), 'AlterationType'] = 'SCNA - AMP'
plotDF[grepl('\\.DEL', rownames(plotDF)), 'AlterationType'] = 'SCNA - DEL'

group.colors <- c("#000000","#E3140F")
alteration.colors = c("#eb004e", "#1cd9d6", "#000000", "#008c59")

p4 <- ggplot(plotDF, aes(x=Train_F, y=Test_F, color=AlterationType, shape=significant_fisher)) + 
  geom_point() +
  xlim(0, 0.45) +
  ylim(0, 0.45) +
  theme_bw() +
  ggtitle(paste('Train/Test Frequencies - Full Features', '(N=', nrow(plotDF) ,')', sep='')) +
  ylab('Test Set Frequency') +
  xlab('Train Set Frequency') +
  theme(axis.text.x = element_text(colour="grey20",size=10,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="grey20",size=10,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="grey20",size=15,angle=0,hjust=.5,vjust=0,face="plain"),
        axis.title.y = element_text(colour="grey20",size=15,hjust=.5,vjust=.5,face="plain")) +
  geom_text(aes(label=ifelse(qval_fisher <= 0.10, as.character(X),'')), vjust=-0.4, size=2) +
  geom_abline() +
  #geom_abline(slope=1, intercept=0.10, linetype='dashed') +
  #geom_abline(slope=1, intercept=-0.10, linetype='dashed') +
  geom_text(x=0.2, y=0.4, 
            label=paste("Fraction significant (qval <= 0.10) = ", 
                        round(sum(plotDF$qval_fisher <= 0.10) / nrow(plotDF), 3), sep=''),
            color='blue') +
  geom_text(x=0.2, y=0.4-0.03, 
            label=paste('pearson correlation:', round(as.numeric(cor.test(plotDF$Train_F, plotDF$Test_F)[4]), 3), '\n',
                        'p =', format(as.numeric(cor.test(plotDF$Train_F, plotDF$Test_F)[3]), scientific=TRUE)),
            color='blue') +
  scale_colour_manual(values = alteration.colors) +
  scale_shape_manual(values = c(19, 4))

plotheight = 6
plotwidth = 7

#ggsave('plots/paper_figures/cohort_frequency_scatter.pdf', p, height = plotheight, width = plotwidth)
#ggsave('plots/cohort_frequency_scatter.png', p, height = plotheight, width = plotwidth)

#ggsave('plots/paper_figures/frequency_vs_coverage.pdf', p_cov, height = plotheight, width = 10)
#ggsave('plots/frequency_vs_coverage.png', p_cov, height = plotheight, width = 10)

ggsave('plots/paper_figures/frequency_plots/cohort_frequency_scatter.pdf', p, height = plotheight, width = plotwidth)
ggsave('plots/frequency_plots/cohort_frequency_scatter.png', p, height = plotheight, width = plotwidth)

ggsave('plots/paper_figures/frequency_plots/cohort_frequency_scatter_natmed.pdf', p1, height = plotheight, width = plotwidth)
ggsave('plots/frequency_plots/cohort_frequency_scatter_natmed.png', p1, height = plotheight, width = plotwidth)

ggsave('plots/paper_figures/cohort_frequency_volcano.pdf', p2, height = plotheight, width = plotwidth)
ggsave('plots/frequency_plots/cohort_frequency_volcano.png', p2, height = plotheight, width = plotwidth)

ggsave('plots/paper_figures/frequency_plots/cohort_frequency_natmed_volcano.pdf', p3, height = plotheight, width = plotwidth)
ggsave('plots/frequency_plots/cohort_frequency_natmed_volcano.png', p3, height = plotheight, width = plotwidth)

ggsave('plots/paper_figures/frequency_plots/traintest_frequency_scatter.pdf', p4, height = plotheight, width = plotwidth)
ggsave('plots/frequency_plots/traintest_frequency_scatter.png', p4, height = plotheight, width = plotwidth)
