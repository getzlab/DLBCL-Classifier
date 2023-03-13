rm(list = ls())

source('src_R/load_libraries.R')

fullDF = read.csv('data_tables/gsm/DLBCL.699.fullGSM.Sep_23_2022.tsv', 
                  sep='\t', row.names = 1)
qvalDF = read.csv('data_tables/qval_dfs/fisher_exact_5x2.Sep_23_2022.combined.tsv', sep='\t', row.names=1)

nmDF = read.csv('data_tables/gsm/DLBCL_significant_event_matrix_NatMed.txt', sep='\t', row.names=1)
rownames(nmDF) = toupper(make.names(rownames(nmDF)))
colnames(nmDF) = toupper(colnames(nmDF))
nmDF = nmDF[ , colSums(is.na(nmDF))==0]
nmlabels = read.csv('data_tables/clustering_labels/NatMed.DLBCL.bestclus.txt',sep='\t', row.names=1)

labels = read.csv('data_tables/confidence_tables/baseline_probabilities.connectivity_based.sensitivity_power2.Sep_23_2022.tsv',
                  sep='\t', row.names=1)

qvalDF = qvalDF[!rownames(qvalDF) %in% c("OR10V1" ,"OR51B6"),]
qvalDF = qvalDF[qvalDF$q <= 0.10,]

cohorts = read.csv('data_tables/sample_sets/sample_inclusion_table.tsv', sep='\t', row.names=1)
cohorts = cohorts[cohorts$included_in_clustering == 'True',]
cohorts = cohorts[rownames(labels), ]

fullDF = data.frame(t(fullDF))
fullDF = fullDF[, !grepl('.CCF', colnames(fullDF))]
fullDF = fullDF[, ]
fullDF = fullDF[, !colnames(fullDF) %in% c('PLOIDY', 'EBV', 'COO')]

for(col in colnames(fullDF)){
  fullDF[, col] = as.numeric(fullDF[, col])
}

NCI_DF = fullDF[rownames(cohorts)[cohorts$cohort == 'Staudt'],]
TCGA_DF = fullDF[rownames(fullDF) %in% rownames(cohorts)[cohorts$cohort == 'TCGA'],]

Staudt_DF = rbind(NCI_DF, TCGA_DF)
Our_DF = fullDF[rownames(fullDF) %in% rownames(cohorts)[cohorts$cohort == 'Shipp'], ]

c1_genes = rownames(qvalDF[((qvalDF$cluster == 'C1') & (qvalDF$q < 0.10)) | ((qvalDF$C1_nf >= 0.30) & (qvalDF$q < 0.10)),])
c2_genes = rownames(qvalDF[((qvalDF$cluster == 'C2') & (qvalDF$q < 0.10)) | ((qvalDF$C2_nf >= 0.30) & (qvalDF$q < 0.10)),])
c3_genes = rownames(qvalDF[((qvalDF$cluster == 'C3') & (qvalDF$q < 0.10)) | ((qvalDF$C3_nf >= 0.30) & (qvalDF$q < 0.10)),])
c4_genes = rownames(qvalDF[((qvalDF$cluster == 'C4') & (qvalDF$q < 0.10)) | ((qvalDF$C4_nf >= 0.30) & (qvalDF$q < 0.10)),])
c5_genes = rownames(qvalDF[((qvalDF$cluster == 'C5') & (qvalDF$q < 0.10)) | ((qvalDF$C5_nf >= 0.30) & (qvalDF$q < 0.10)),])

nm_c1 = rownames(nmDF)[rownames(nmDF) %in% c1_genes]
nm_c2 = rownames(nmDF)[rownames(nmDF) %in% c2_genes]
nm_c3 = rownames(nmDF)[rownames(nmDF) %in% c3_genes]
nm_c4 = rownames(nmDF)[rownames(nmDF) %in% c4_genes]
nm_c5 = rownames(nmDF)[rownames(nmDF) %in% c5_genes]

c1_x = intersect(nm_c1, c1_genes)
c2_x = intersect(nm_c2, c2_genes)
c3_x = intersect(nm_c3, c3_genes)
c4_x = intersect(nm_c4, c4_genes)
c5_x = intersect(nm_c5, c5_genes)

c1_st_s = rownames(labels[(labels$cluster == 1) & (rownames(labels) %in% rownames(Staudt_DF)),])
c1_sh_s = rownames(labels[(labels$cluster == 1) & (rownames(labels) %in% rownames(Our_DF)),])
c2_st_s = rownames(labels[(labels$cluster == 2) & (rownames(labels) %in% rownames(Staudt_DF)),])
c2_sh_s = rownames(labels[(labels$cluster == 2) & (rownames(labels) %in% rownames(Our_DF)),])
c3_st_s = rownames(labels[(labels$cluster == 3) & (rownames(labels) %in% rownames(Staudt_DF)),])
c3_sh_s = rownames(labels[(labels$cluster == 3) & (rownames(labels) %in% rownames(Our_DF)),])
c4_st_s = rownames(labels[(labels$cluster == 4) & (rownames(labels) %in% rownames(Staudt_DF)),])
c4_sh_s = rownames(labels[(labels$cluster == 4) & (rownames(labels) %in% rownames(Our_DF)),])
c5_st_s = rownames(labels[(labels$cluster == 5) & (rownames(labels) %in% rownames(Staudt_DF)),])
c5_sh_s = rownames(labels[(labels$cluster == 5) & (rownames(labels) %in% rownames(Our_DF)),])

c1_nm_s = rownames(nmlabels[nmlabels$cluster == 1, ])
c2_nm_s = rownames(nmlabels[nmlabels$cluster == 2, ])
c3_nm_s = rownames(nmlabels[nmlabels$cluster == 3, ])
c4_nm_s = rownames(nmlabels[nmlabels$cluster == 4, ])
c5_nm_s = rownames(nmlabels[nmlabels$cluster == 5, ])

c1_nm_s = c1_nm_s[c1_nm_s %in% colnames(nmDF)]
c2_nm_s = c2_nm_s[c2_nm_s %in% colnames(nmDF)]
c3_nm_s = c3_nm_s[c3_nm_s %in% colnames(nmDF)]
c4_nm_s = c4_nm_s[c4_nm_s %in% colnames(nmDF)]
c5_nm_s = c5_nm_s[c5_nm_s %in% colnames(nmDF)]

nmDF = data.frame(t(nmDF))


plotDF_c1 = data.frame(cbind(colSums(Staudt_DF[c1_st_s, c1_x] != 0) / nrow(Staudt_DF[c1_st_s, c1_x]), 
                             colSums(Our_DF[c1_sh_s, c1_x] != 0) / nrow(Our_DF[c1_sh_s, c1_x]),
                             colSums(nmDF[c1_nm_s, c1_x] != 0) / nrow(nmDF[c1_nm_s, c1_x]),
                             colSums(fullDF[c(c1_st_s, c1_sh_s), c1_x] != 0) / nrow(fullDF[c(c1_st_s, c1_sh_s), c1_x]),
                             colSums(Staudt_DF[c1_st_s, c1_x] != 0), 
                             colSums(Our_DF[c1_sh_s, c1_x] != 0),
                             colSums(nmDF[c1_nm_s, c1_x] != 0),
                             colSums(fullDF[c(c1_st_s, c1_sh_s), c1_x] != 0),
                             nrow(Staudt_DF[c1_st_s, c1_x]),
                             nrow(Our_DF[c1_sh_s, c1_x]),
                             nrow(nmDF[c1_nm_s, c1_x]),
                             nrow(fullDF[c(c1_st_s, c1_sh_s), c1_x])))

plotDF_c2 = data.frame(cbind(colSums(Staudt_DF[c2_st_s, c2_x] != 0) / nrow(Staudt_DF[c2_st_s, c2_x]), 
                             colSums(Our_DF[c2_sh_s, c2_x] != 0) / nrow(Our_DF[c2_sh_s, c2_x]),
                             colSums(nmDF[c2_nm_s, c2_x] != 0) / nrow(nmDF[c2_nm_s, c2_x]),
                             colSums(fullDF[c(c2_st_s, c2_sh_s), c2_x] != 0) / nrow(fullDF[c(c2_st_s, c2_sh_s), c2_x]),
                             colSums(Staudt_DF[c2_st_s, c2_x] != 0), 
                             colSums(Our_DF[c2_sh_s, c2_x] != 0),
                             colSums(nmDF[c2_nm_s, c2_x] != 0),
                             colSums(fullDF[c(c2_st_s, c2_sh_s), c2_x] != 0),
                             nrow(Staudt_DF[c2_st_s, c2_x]),
                             nrow(Our_DF[c2_sh_s, c2_x]),
                             nrow(nmDF[c2_nm_s, c2_x]),
                             nrow(fullDF[c(c2_st_s, c2_sh_s), c2_x])))

plotDF_c3 = data.frame(cbind(colSums(Staudt_DF[c3_st_s, c3_x] != 0) / nrow(Staudt_DF[c3_st_s, c3_x]), 
                             colSums(Our_DF[c3_sh_s, c3_x] != 0) / nrow(Our_DF[c3_sh_s, c3_x]),
                             colSums(nmDF[c3_nm_s, c3_x] != 0) / nrow(nmDF[c3_nm_s, c3_x]),
                             colSums(fullDF[c(c3_st_s, c3_sh_s), c3_x] != 0) / nrow(fullDF[c(c3_st_s, c3_sh_s), c3_x]),
                             colSums(Staudt_DF[c3_st_s, c3_x] != 0), 
                             colSums(Our_DF[c3_sh_s, c3_x] != 0),
                             colSums(nmDF[c3_nm_s, c3_x] != 0),
                             colSums(fullDF[c(c3_st_s, c3_sh_s), c3_x] != 0),
                             nrow(Staudt_DF[c3_st_s, c3_x]),
                             nrow(Our_DF[c3_sh_s, c3_x]),
                             nrow(nmDF[c3_nm_s, c3_x]),
                             nrow(fullDF[c(c3_st_s, c3_sh_s), c3_x])))

plotDF_c4 = data.frame(data.frame(cbind(colSums(Staudt_DF[c4_st_s, c4_x] != 0) / nrow(Staudt_DF[c4_st_s, c4_x]), 
                                        colSums(Our_DF[c4_sh_s, c4_x] != 0) / nrow(Our_DF[c4_sh_s, c4_x]),
                                        colSums(nmDF[c4_nm_s, c4_x] != 0) / nrow(nmDF[c4_nm_s, c4_x]),
                                        colSums(fullDF[c(c4_st_s, c4_sh_s), c4_x] != 0) / nrow(fullDF[c(c4_st_s, c4_sh_s), c4_x]),
                                        colSums(Staudt_DF[c4_st_s, c4_x] != 0), 
                                        colSums(Our_DF[c4_sh_s, c4_x] != 0),
                                        colSums(nmDF[c4_nm_s, c4_x] != 0),
                                        colSums(fullDF[c(c4_st_s, c4_sh_s), c4_x] != 0),
                                        nrow(Staudt_DF[c4_st_s, c4_x]),
                                        nrow(Our_DF[c4_sh_s, c4_x]),
                                        nrow(nmDF[c4_nm_s, c4_x]),
                                        nrow(fullDF[c(c4_st_s, c4_sh_s), c4_x]))))

plotDF_c5 = data.frame(data.frame(cbind(colSums(Staudt_DF[c5_st_s, c5_x] != 0) / nrow(Staudt_DF[c5_st_s, c5_x]), 
                                        colSums(Our_DF[c5_sh_s, c5_x] != 0) / nrow(Our_DF[c5_sh_s, c5_x]),
                                        colSums(nmDF[c5_nm_s, c5_x] != 0) / nrow(nmDF[c5_nm_s, c5_x]),
                                        colSums(fullDF[c(c5_st_s, c5_sh_s), c5_x] != 0) / nrow(fullDF[c(c5_st_s, c5_sh_s), c5_x]),
                                        colSums(Staudt_DF[c5_st_s, c5_x] != 0), 
                                        colSums(Our_DF[c5_sh_s, c5_x] != 0),
                                        colSums(nmDF[c5_nm_s, c5_x] != 0),
                                        colSums(fullDF[c(c5_st_s, c5_sh_s), c5_x] != 0),
                                        nrow(Staudt_DF[c5_st_s, c5_x]),
                                        nrow(Our_DF[c5_sh_s, c5_x]),
                                        nrow(nmDF[c5_nm_s, c5_x]),
                                        nrow(fullDF[c(c5_st_s, c5_sh_s), c5_x]))))

colnames(plotDF_c1) = c('Schmitz_F', 'Chapuy_F', 'NM_F', 'Combined_F',
                        'Schmitz_c', 'Chapuy_c', 'NM_c', 'Combined_c',
                        'Schmitz_n', 'Chapuy_n', 'NM_n', 'Combined_n')

plotDF_c1$Schmitz_lower_sig = apply(plotDF_c1, 1, 
                            function(x) binom.test(as.numeric(x[5]), as.numeric(x[9]), as.numeric(x[1]), conf.level = 0.68)[[4]][[1]])
plotDF_c1$Schmitz_upper_sig = apply(plotDF_c1, 1, 
                            function(x) binom.test(as.numeric(x[5]), as.numeric(x[9]), as.numeric(x[1]), conf.level = 0.68)[[4]][[2]])

colnames(plotDF_c2) = c('Schmitz_F', 'Chapuy_F', 'NM_F', 'Combined_F',
                        'Schmitz_c', 'Chapuy_c', 'NM_c', 'Combined_c',
                        'Schmitz_n', 'Chapuy_n', 'NM_n', 'Combined_n')

colnames(plotDF_c3) = c('Schmitz_F', 'Chapuy_F', 'NM_F', 'Combined_F',
                        'Schmitz_c', 'Chapuy_c', 'NM_c', 'Combined_c',
                        'Schmitz_n', 'Chapuy_n', 'NM_n', 'Combined_n')

colnames(plotDF_c4) = c('Schmitz_F', 'Chapuy_F', 'NM_F', 'Combined_F',
                        'Schmitz_c', 'Chapuy_c', 'NM_c', 'Combined_c',
                        'Schmitz_n', 'Chapuy_n', 'NM_n', 'Combined_n')

colnames(plotDF_c5) = c('Schmitz_F', 'Chapuy_F', 'NM_F', 'Combined_F',
                        'Schmitz_c', 'Chapuy_c', 'NM_c', 'Combined_c',
                        'Schmitz_n', 'Chapuy_n', 'NM_n', 'Combined_n')

plotDF_c1$AlterationType = 'SNV'
plotDF_c1[grepl('SV', rownames(plotDF_c1)), 'AlterationType'] = 'SV'
plotDF_c1[grepl('\\.AMP', rownames(plotDF_c1)), 'AlterationType'] = 'SCNA - AMP'
plotDF_c1[grepl('\\.DEL', rownames(plotDF_c1)), 'AlterationType'] = 'SCNA - DEL'

plotDF_c2$AlterationType = 'SNV'
plotDF_c2[grepl('SV', rownames(plotDF_c2)), 'AlterationType'] = 'SV'
plotDF_c2[grepl('\\.AMP', rownames(plotDF_c2)), 'AlterationType'] = 'SCNA - AMP'
plotDF_c2[grepl('\\.DEL', rownames(plotDF_c2)), 'AlterationType'] = 'SCNA - DEL'

plotDF_c3$AlterationType = 'SNV'
plotDF_c3[grepl('SV', rownames(plotDF_c3)), 'AlterationType'] = 'SV'
plotDF_c3[grepl('\\.AMP', rownames(plotDF_c3)), 'AlterationType'] = 'SCNA - AMP'
plotDF_c3[grepl('\\.DEL', rownames(plotDF_c3)), 'AlterationType'] = 'SCNA - DEL'

plotDF_c4$AlterationType = 'SNV'
plotDF_c4[grepl('SV', rownames(plotDF_c4)), 'AlterationType'] = 'SV'
plotDF_c4[grepl('\\.AMP', rownames(plotDF_c4)), 'AlterationType'] = 'SCNA - AMP'
plotDF_c4[grepl('\\.DEL', rownames(plotDF_c4)), 'AlterationType'] = 'SCNA - DEL'

plotDF_c5$AlterationType = 'SNV'
plotDF_c5[grepl('SV', rownames(plotDF_c5)), 'AlterationType'] = 'SV'
plotDF_c5[grepl('\\.AMP', rownames(plotDF_c5)), 'AlterationType'] = 'SCNA - AMP'
plotDF_c5[grepl('\\.DEL', rownames(plotDF_c5)), 'AlterationType'] = 'SCNA - DEL'

alteration_colors_c1 = c("#eb004e", "#1cd9d6", "#000000", "#008c59")
alteration_colors_c2 = c("#eb004e", "#1cd9d6", "#000000")
alteration_colors_c3 = c("#1cd9d6", "#000000", "#008c59")
alteration_colors_c4 = c("#000000")
alteration_colors_c5 = c("#eb004e", "#1cd9d6", "#000000")


point_size = 3
alpha = 0.65

##########################
# C1 -> C5 nm vs schmitz #
##########################

p_c1_st <- ggplot(plotDF_c1, aes(x=NM_F, y=Schmitz_F, color=AlterationType)) + 
  geom_point(size=point_size, alpha=alpha) +
  xlim(c(0, 0.8)) +
  ylim(c(0, 0.8)) +
  theme_bw() +
  ggtitle(paste('NM vs Schmitz', '(N=', nrow(plotDF_c1) ,')', sep='')) +
  ylab('Schmitz et al. Frequencies') +
  xlab('Nature Med Frequencies') +
  theme(axis.text.x = element_text(colour="grey20",size=10,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="grey20",size=10,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="grey20",size=15,angle=0,hjust=.5,vjust=0,face="plain"),
        axis.title.y = element_text(colour="grey20",size=15,hjust=.5,vjust=.5,face="plain")) +
  #geom_text(aes(label=rownames(plotDF_c1),vjust=-0.4), show.legend = F) +
  geom_abline() +
  geom_text(x=0.2, y = 0.7,
            label=paste('Pearson correlation:', round(as.numeric(cor.test(plotDF_c1$NM_F, plotDF_c1$Schmitz_F)[4]), 2)), color='blue') +
  scale_colour_manual(values = alteration_colors_c1)

p_c2_st <- ggplot(plotDF_c2, aes(x=NM_F, y=Schmitz_F, color=AlterationType)) + 
  geom_point(size=point_size, alpha=alpha) +
  xlim(c(0, 0.8)) +
  ylim(c(0, 0.8)) +
  theme_bw() +
  ggtitle(paste('NM vs Schmitz', '(N=', nrow(plotDF_c2) ,')', sep='')) +
  ylab('Schmitz et al. Frequencies') +
  xlab('Nature Med Frequencies') +
  theme(axis.text.x = element_text(colour="grey20",size=10,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="grey20",size=10,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="grey20",size=15,angle=0,hjust=.5,vjust=0,face="plain"),
        axis.title.y = element_text(colour="grey20",size=15,hjust=.5,vjust=.5,face="plain")) +
  #geom_text(aes(label=rownames(plotDF_c2),vjust=-0.4, size=2)) +
  geom_abline() +
  geom_text(x=0.2, y = 0.7,
            label=paste('Pearson correlation:', round(as.numeric(cor.test(plotDF_c2$NM_F, plotDF_c2$Schmitz_F)[4]), 2)),
            color='blue') +
  scale_colour_manual(values = alteration_colors_c2)

p_c3_st <- ggplot(plotDF_c3, aes(x=NM_F, y=Schmitz_F, color=AlterationType)) + 
  geom_point(size=point_size, alpha=alpha) +
  xlim(c(0, 0.8)) +
  ylim(c(0, 0.8)) +
  theme_bw() +
  ggtitle(paste('NM vs Schmitz', '(N=', nrow(plotDF_c3) ,')', sep='')) +
  ylab('Schmitz et al. Frequencies') +
  xlab('Nature Med Frequencies') +
  theme(axis.text.x = element_text(colour="grey20",size=10,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="grey20",size=10,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="grey20",size=15,angle=0,hjust=.5,vjust=0,face="plain"),
        axis.title.y = element_text(colour="grey20",size=15,hjust=.5,vjust=.5,face="plain")) +
  #geom_text(aes(label=rownames(plotDF_c3),vjust=-0.4, size=2)) +
  geom_abline() +
  geom_text(x=0.2, y = 0.7,
            label=paste('Pearson correlation:', round(as.numeric(cor.test(plotDF_c3$NM_F, plotDF_c3$Schmitz_F)[4]), 2)),
            color='blue') +
  scale_colour_manual(values = alteration_colors_c3)

p_c4_st <- ggplot(plotDF_c4, aes(x=NM_F, y=Schmitz_F, color=AlterationType)) + 
  geom_point(size=point_size, alpha=alpha) +
  xlim(c(0, 0.8)) +
  ylim(c(0, 0.8)) +
  theme_bw() +
  ggtitle(paste('NM vs Schmitz', '(N=', nrow(plotDF_c4) ,')', sep='')) +
  ylab('Schmitz et al. Frequencies') +
  xlab('Nature Med Frequencies') +
  theme(axis.text.x = element_text(colour="grey20",size=10,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="grey20",size=10,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="grey20",size=15,angle=0,hjust=.5,vjust=0,face="plain"),
        axis.title.y = element_text(colour="grey20",size=15,hjust=.5,vjust=.5,face="plain")) +
  #geom_text(aes(label=rownames(plotDF_c4),vjust=-0.4, size=2)) +
  geom_abline() +
  geom_text(x=0.2, y = 0.7,
            label=paste('Pearson correlation:', round(as.numeric(cor.test(plotDF_c4$NM_F, plotDF_c4$Schmitz_F)[4]), 2)),
            color='blue') +
  scale_colour_manual(values = alteration_colors_c4)

p_c5_st <- ggplot(plotDF_c5, aes(x=NM_F, y=Schmitz_F, color=AlterationType)) + 
  geom_point(size=point_size, alpha=alpha) +
  xlim(c(0, 0.8)) +
  ylim(c(0, 0.8)) +
  theme_bw() +
  ggtitle(paste('NM vs Schmitz', '(N=', nrow(plotDF_c5) ,')', sep='')) +
  ylab('Schmitz et al. Frequencies') +
  xlab('Nature Med Frequencies') +
  theme(axis.text.x = element_text(colour="grey20",size=10,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="grey20",size=10,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="grey20",size=15,angle=0,hjust=.5,vjust=0,face="plain"),
        axis.title.y = element_text(colour="grey20",size=15,hjust=.5,vjust=.5,face="plain")) +
  #geom_text(aes(label=rownames(plotDF_c5),vjust=-0.4, size=2)) +
  geom_abline() +
  geom_text(x=0.2, y = 0.7,
            label=paste('Pearson correlation:', round(as.numeric(cor.test(plotDF_c5$NM_F, plotDF_c5$Schmitz_F)[4]),2)),
            color='blue') +
  scale_colour_manual(values = alteration_colors_c5)

##########################
# C1 -> C5 nm vs ours    #
##########################

p_c1_sh <- ggplot(plotDF_c1, aes(x=NM_F, y=Chapuy_F, color=AlterationType)) + 
  geom_point(size=point_size, alpha=alpha) +
  xlim(c(0, 0.8)) +
  ylim(c(0, 0.8)) +
  theme_bw() +
  ggtitle(paste('NM vs Ours', '(N=', nrow(plotDF_c1) ,')', sep='')) +
  ylab('Our Frequencies') +
  xlab('Nature Med Frequencies') +
  theme(axis.text.x = element_text(colour="grey20",size=10,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="grey20",size=10,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="grey20",size=15,angle=0,hjust=.5,vjust=0,face="plain"),
        axis.title.y = element_text(colour="grey20",size=15,hjust=.5,vjust=.5,face="plain")) +
  #geom_text(aes(label=rownames(plotDF_c1),vjust=-0.4, size=2)) +
  geom_abline() +
  geom_text(x=0.2, y = 0.7,
            label=paste('Pearson correlation:', round(as.numeric(cor.test(plotDF_c1$NM_F, plotDF_c1$Chapuy_F)[4]), 2)),
            color='blue') +
  scale_colour_manual(values = alteration_colors_c1)

p_c2_sh <- ggplot(plotDF_c2, aes(x=NM_F, y=Chapuy_F, color=AlterationType)) + 
  geom_point(size=point_size, alpha=alpha) +
  xlim(c(0, 0.8)) +
  ylim(c(0, 0.8)) +
  theme_bw() +
  ggtitle(paste('NM vs Ours', '(N=', nrow(plotDF_c2) ,')', sep='')) +
  ylab('Our Frequencies') +
  xlab('Nature Med Frequencies') +
  theme(axis.text.x = element_text(colour="grey20",size=10,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="grey20",size=10,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="grey20",size=15,angle=0,hjust=.5,vjust=0,face="plain"),
        axis.title.y = element_text(colour="grey20",size=15,hjust=.5,vjust=.5,face="plain")) +
  #geom_text(aes(label=rownames(plotDF_c2),vjust=-0.4, size=2)) +
  geom_abline() +
  geom_text(x=0.2, y = 0.7,
            label=paste('Pearson correlation:', round(as.numeric(cor.test(plotDF_c2$NM_F, plotDF_c2$Chapuy_F)[4]), 2)),
            color='blue') +
  scale_colour_manual(values = alteration_colors_c2)

p_c3_sh <- ggplot(plotDF_c3, aes(x=NM_F, y=Chapuy_F, color=AlterationType)) + 
  geom_point(size=point_size, alpha=alpha) +
  xlim(c(0, 0.8)) +
  ylim(c(0, 0.8)) +
  theme_bw() +
  ggtitle(paste('NM vs Ours', '(N=', nrow(plotDF_c3) ,')', sep='')) +
  ylab('Our Frequencies') +
  xlab('Nature Med Frequencies') +
  theme(axis.text.x = element_text(colour="grey20",size=10,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="grey20",size=10,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="grey20",size=15,angle=0,hjust=.5,vjust=0,face="plain"),
        axis.title.y = element_text(colour="grey20",size=15,hjust=.5,vjust=.5,face="plain")) +
  #geom_text(aes(label=rownames(plotDF_c3),vjust=-0.4, size=2)) +
  geom_abline() +
  geom_text(x=0.2, y = 0.7,
            label=paste('Pearson correlation:', round(as.numeric(cor.test(plotDF_c3$NM_F, plotDF_c3$Chapuy_F)[4]), 2)),
            color='blue') +
  scale_colour_manual(values = alteration_colors_c3)

p_c4_sh <- ggplot(plotDF_c4, aes(x=NM_F, y=Chapuy_F, color=AlterationType)) + 
  geom_point(size=point_size, alpha=alpha) +
  xlim(c(0, 0.8)) +
  ylim(c(0, 0.8)) +
  theme_bw() +
  ggtitle(paste('NM vs Ours', '(N=', nrow(plotDF_c4) ,')', sep='')) +
  ylab('Our Frequencies') +
  xlab('Nature Med Frequencies') +
  theme(axis.text.x = element_text(colour="grey20",size=10,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="grey20",size=10,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="grey20",size=15,angle=0,hjust=.5,vjust=0,face="plain"),
        axis.title.y = element_text(colour="grey20",size=15,hjust=.5,vjust=.5,face="plain")) +
  #geom_text(aes(label=rownames(plotDF_c4),vjust=-0.4, size=2)) +
  geom_abline() +
  geom_text(x=0.2, y = 0.7,
            label=paste('Pearson correlation:', round(as.numeric(cor.test(plotDF_c4$NM_F, plotDF_c4$Chapuy_F)[4]), 2)),
            color='blue') +
  scale_colour_manual(values = alteration_colors_c4)

p_c5_sh <- ggplot(plotDF_c5, aes(x=NM_F, y=Chapuy_F, color=AlterationType)) + 
  geom_point(size=point_size, alpha=alpha) +
  xlim(c(0, 0.8)) +
  ylim(c(0, 0.8)) +
  theme_bw() +
  ggtitle(paste('NM vs Ours', '(N=', nrow(plotDF_c5) ,')', sep='')) +
  ylab('Our Frequencies') +
  xlab('Nature Med Frequencies') +
  theme(axis.text.x = element_text(colour="grey20",size=10,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="grey20",size=10,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="grey20",size=15,angle=0,hjust=.5,vjust=0,face="plain"),
        axis.title.y = element_text(colour="grey20",size=15,hjust=.5,vjust=.5,face="plain")) +
  #geom_text(aes(label=rownames(plotDF_c5),vjust=-0.4, size=2)) +
  geom_abline() +
  geom_text(x=0.2, y = 0.7,
            label=paste('Pearson correlation:', round(as.numeric(cor.test(plotDF_c5$NM_F, plotDF_c5$Chapuy_F)[4]), 2)),
            color='blue') +
  scale_colour_manual(values = alteration_colors_c5)

##########################
# C1 -> C5 nm vs comb    #
##########################

p_c1_c <- ggplot(plotDF_c1, aes(x=NM_F, y=Combined_F, color=AlterationType)) + 
  geom_point(size=point_size, alpha=alpha) +
  xlim(c(0, 0.8)) +
  ylim(c(0, 0.8)) +
  theme_bw() +
  ggtitle(paste('NM vs Combined', '(N=', nrow(plotDF_c1) ,')', sep='')) +
  ylab('Combined Frequency') +
  xlab('Nature Med Frequencies') +
  theme(axis.text.x = element_text(colour="grey20",size=10,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="grey20",size=10,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="grey20",size=15,angle=0,hjust=.5,vjust=0,face="plain"),
        axis.title.y = element_text(colour="grey20",size=15,hjust=.5,vjust=.5,face="plain")) +
  #geom_text(aes(label=rownames(plotDF_c1),vjust=-0.4, size=2)) +
  geom_abline() +
  geom_text(x=0.2, y = 0.7,
            label=paste('Pearson correlation:', round(as.numeric(cor.test(plotDF_c1$NM_F, plotDF_c1$Combined_F)[4]), 2)),
            color='blue') +
  scale_colour_manual(values = alteration_colors_c1)

p_c2_c <- ggplot(plotDF_c2, aes(x=NM_F, y=Combined_F, color=AlterationType)) + 
  geom_point(size=point_size, alpha=alpha) +
  xlim(c(0, 0.8)) +
  ylim(c(0, 0.8)) +
  theme_bw() +
  ggtitle(paste('NM vs Combined', '(N=', nrow(plotDF_c2) ,')', sep='')) +
  ylab('Combined Frequency') +
  xlab('Nature Med Frequencies') +
  theme(axis.text.x = element_text(colour="grey20",size=10,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="grey20",size=10,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="grey20",size=15,angle=0,hjust=.5,vjust=0,face="plain"),
        axis.title.y = element_text(colour="grey20",size=15,hjust=.5,vjust=.5,face="plain")) +
  #geom_text(aes(label=rownames(plotDF_c2),vjust=-0.4, size=2)) +
  geom_abline() +
  geom_text(x=0.2, y = 0.7,
            label=paste('Pearson correlation:', round(as.numeric(cor.test(plotDF_c2$NM_F, plotDF_c2$Combined_F)[4]), 2)),
            color='blue') +
  scale_colour_manual(values = alteration_colors_c2)

p_c3_c <- ggplot(plotDF_c3, aes(x=NM_F, y=Combined_F, color=AlterationType)) + 
  geom_point(size=point_size, alpha=alpha) +
  xlim(c(0, 0.8)) +
  ylim(c(0, 0.8)) +
  theme_bw() +
  ggtitle(paste('NM vs Combined', '(N=', nrow(plotDF_c3) ,')', sep='')) +
  ylab('Combined Frequency') +
  xlab('Nature Med Frequencies') +
  theme(axis.text.x = element_text(colour="grey20",size=10,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="grey20",size=10,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="grey20",size=15,angle=0,hjust=.5,vjust=0,face="plain"),
        axis.title.y = element_text(colour="grey20",size=15,hjust=.5,vjust=.5,face="plain")) +
  #geom_text(aes(label=rownames(plotDF_c3),vjust=-0.4, size=2)) +
  geom_abline() +
  geom_text(x=0.2, y = 0.7,
            label=paste('Pearson correlation:', round(as.numeric(cor.test(plotDF_c3$NM_F, plotDF_c3$Combined_F)[4]), 2)),
            color='blue') +
  scale_colour_manual(values = alteration_colors_c3)

p_c4_c <- ggplot(plotDF_c4, aes(x=NM_F, y=Combined_F, color=AlterationType)) + 
  geom_point(size=point_size, alpha=alpha) +
  xlim(c(0, 0.8)) +
  ylim(c(0, 0.8)) +
  theme_bw() +
  ggtitle(paste('NM vs Combined', '(N=', nrow(plotDF_c4) ,')', sep='')) +
  ylab('Combined Frequency') +
  xlab('Nature Med Frequencies') +
  theme(axis.text.x = element_text(colour="grey20",size=10,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="grey20",size=10,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="grey20",size=15,angle=0,hjust=.5,vjust=0,face="plain"),
        axis.title.y = element_text(colour="grey20",size=15,hjust=.5,vjust=.5,face="plain")) +
  #geom_text(aes(label=rownames(plotDF_c4),vjust=-0.4, size=2)) +
  geom_abline() +
  geom_text(x=0.2, y = 0.7,
            label=paste('Pearson correlation:', round(as.numeric(cor.test(plotDF_c4$NM_F, plotDF_c4$Combined_F)[4]), 2)),
            color='blue') +
  scale_colour_manual(values = alteration_colors_c4)

p_c5_c <- ggplot(plotDF_c5, aes(x=NM_F, y=Combined_F, color=AlterationType)) + 
  geom_point(size=point_size, alpha=alpha) +
  xlim(c(0, 0.8)) +
  ylim(c(0, 0.8)) +
  theme_bw() +
  ggtitle(paste('NM vs Combined', '(N=', nrow(plotDF_c5) ,')', sep='')) +
  ylab('Combined Frequency') +
  xlab('Nature Med Frequencies') +
  theme(axis.text.x = element_text(colour="grey20",size=10,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="grey20",size=10,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="grey20",size=15,angle=0,hjust=.5,vjust=0,face="plain"),
        axis.title.y = element_text(colour="grey20",size=15,hjust=.5,vjust=.5,face="plain")) +
  #geom_text(aes(label=rownames(plotDF_c5),vjust=-0.4, size=2)) +
  geom_abline() +
  geom_text(x=0.2, y = 0.7,
            label=paste('Pearson correlation:', round(as.numeric(cor.test(plotDF_c5$NM_F, plotDF_c5$Combined_F)[4]),2)),
            color='blue') +
  scale_colour_manual(values = alteration_colors_c5)

plotheight = 6
plotwidth = 7

#####################
# save staudt plots #
#####################

ggsave('plots/paper_figures/frequency_plots/cluster_frequencies/c1_st_nm_nolabels.pdf', p_c1_st, height = plotheight, width = plotwidth)
ggsave('plots/frequency_plots/cluster_frequencies/c1_st_nm_nolabels.png', p_c1_st, height = plotheight, width = plotwidth)
ggsave('plots/paper_figures/frequency_plots/cluster_frequencies/c2_st_nm_nolabels.pdf', p_c2_st, height = plotheight, width = plotwidth)
ggsave('plots/frequency_plots/cluster_frequencies/c2_st_nm_nolabels.png', p_c2_st, height = plotheight, width = plotwidth)
ggsave('plots/paper_figures/frequency_plots/cluster_frequencies/c3_st_nm_nolabels.pdf', p_c3_st, height = plotheight, width = plotwidth)
ggsave('plots/frequency_plots/cluster_frequencies/c3_st_nm_nolabels.png', p_c3_st, height = plotheight, width = plotwidth)
ggsave('plots/paper_figures/frequency_plots/cluster_frequencies/c4_st_nm_nolabels.pdf', p_c4_st, height = plotheight, width = plotwidth)
ggsave('plots/frequency_plots/cluster_frequencies/c4_st_nm_nolabels.png', p_c4_st, height = plotheight, width = plotwidth)
ggsave('plots/paper_figures/frequency_plots/cluster_frequencies/c5_st_nm_nolabels.pdf', p_c5_st, height = plotheight, width = plotwidth)
ggsave('plots/frequency_plots/cluster_frequencies/c5_st_nm_nolabels.png', p_c5_st, height = plotheight, width = plotwidth)

p_c1_st = p_c1_st + geom_text(aes(label=rownames(plotDF_c1),vjust=-0.4, size=2), show.legend = F)
p_c2_st = p_c2_st + geom_text(aes(label=rownames(plotDF_c2),vjust=-0.4, size=2), show.legend = F)
p_c3_st = p_c3_st + geom_text(aes(label=rownames(plotDF_c3),vjust=-0.4, size=2), show.legend = F)
p_c4_st = p_c4_st + geom_text(aes(label=rownames(plotDF_c4),vjust=-0.4, size=2), show.legend = F)
p_c5_st = p_c5_st + geom_text(aes(label=rownames(plotDF_c5),vjust=-0.4, size=2), show.legend = F)

ggsave('plots/paper_figures/frequency_plots/cluster_frequencies/c1_st_nm.pdf', p_c1_st, height = plotheight, width = plotwidth)
ggsave('plots/frequency_plots/cluster_frequencies/c1_st_nm.png', p_c1_st, height = plotheight, width = plotwidth)
ggsave('plots/paper_figures/frequency_plots/cluster_frequencies/c2_st_nm.pdf', p_c2_st, height = plotheight, width = plotwidth)
ggsave('plots/frequency_plots/cluster_frequencies/c2_st_nm.png', p_c2_st, height = plotheight, width = plotwidth)
ggsave('plots/paper_figures/frequency_plots/cluster_frequencies/c3_st_nm.pdf', p_c3_st, height = plotheight, width = plotwidth)
ggsave('plots/frequency_plots/cluster_frequencies/c3_st_nm.png', p_c3_st, height = plotheight, width = plotwidth)
ggsave('plots/paper_figures/frequency_plots/cluster_frequencies/c4_st_nm.pdf', p_c4_st, height = plotheight, width = plotwidth)
ggsave('plots/frequency_plots/cluster_frequencies/c4_st_nm.png', p_c4_st, height = plotheight, width = plotwidth)
ggsave('plots/paper_figures/frequency_plots/cluster_frequencies/c5_st_nm.pdf', p_c5_st, height = plotheight, width = plotwidth)
ggsave('plots/frequency_plots/cluster_frequencies/c5_st_nm.png', p_c5_st, height = plotheight, width = plotwidth)

####################
# save shipp plots #
####################

ggsave('plots/paper_figures/frequency_plots/cluster_frequencies/c1_sh_nm_nolabels.pdf', p_c1_sh, height = plotheight, width = plotwidth)
ggsave('plots/frequency_plots/cluster_frequencies/c1_sh_nm_nolabels.png', p_c1_sh, height = plotheight, width = plotwidth)
ggsave('plots/paper_figures/frequency_plots/cluster_frequencies/c2_sh_nm_nolabels.pdf', p_c2_sh, height = plotheight, width = plotwidth)
ggsave('plots/frequency_plots/cluster_frequencies/c2_sh_nm_nolabels.png', p_c2_sh, height = plotheight, width = plotwidth)
ggsave('plots/paper_figures/frequency_plots/cluster_frequencies/c3_sh_nm_nolabels.pdf', p_c3_sh, height = plotheight, width = plotwidth)
ggsave('plots/frequency_plots/cluster_frequencies/c3_sh_nm_nolabels.png', p_c3_sh, height = plotheight, width = plotwidth)
ggsave('plots/paper_figures/frequency_plots/cluster_frequencies/c4_sh_nm_nolabels.pdf', p_c4_sh, height = plotheight, width = plotwidth)
ggsave('plots/frequency_plots/cluster_frequencies/c4_sh_nm_nolabels.png', p_c4_sh, height = plotheight, width = plotwidth)
ggsave('plots/paper_figures/frequency_plots/cluster_frequencies/c5_sh_nm_nolabels.pdf', p_c5_sh, height = plotheight, width = plotwidth)
ggsave('plots/frequency_plots/cluster_frequencies/c5_sh_nm_nolabels.png', p_c5_sh, height = plotheight, width = plotwidth)

p_c1_sh = p_c1_sh + geom_text(aes(label=rownames(plotDF_c1),vjust=-0.4, size=2), show.legend = F)
p_c2_sh = p_c2_sh + geom_text(aes(label=rownames(plotDF_c2),vjust=-0.4, size=2), show.legend = F)
p_c3_sh = p_c3_sh + geom_text(aes(label=rownames(plotDF_c3),vjust=-0.4, size=2), show.legend = F)
p_c4_sh = p_c4_sh + geom_text(aes(label=rownames(plotDF_c4),vjust=-0.4, size=2), show.legend = F)
p_c5_sh = p_c5_sh + geom_text(aes(label=rownames(plotDF_c5),vjust=-0.4, size=2), show.legend = F)

ggsave('plots/paper_figures/frequency_plots/cluster_frequencies/c1_sh_nm.pdf', p_c1_sh, height = plotheight, width = plotwidth)
ggsave('plots/frequency_plots/cluster_frequencies/c1_sh_nm.png', p_c1_sh, height = plotheight, width = plotwidth)
ggsave('plots/paper_figures/frequency_plots/cluster_frequencies/c2_sh_nm.pdf', p_c2_sh, height = plotheight, width = plotwidth)
ggsave('plots/frequency_plots/cluster_frequencies/c2_sh_nm.png', p_c2_sh, height = plotheight, width = plotwidth)
ggsave('plots/paper_figures/frequency_plots/cluster_frequencies/c3_sh_nm.pdf', p_c3_sh, height = plotheight, width = plotwidth)
ggsave('plots/frequency_plots/cluster_frequencies/c3_sh_nm.png', p_c3_sh, height = plotheight, width = plotwidth)
ggsave('plots/paper_figures/frequency_plots/cluster_frequencies/c4_sh_nm.pdf', p_c4_sh, height = plotheight, width = plotwidth)
ggsave('plots/frequency_plots/cluster_frequencies/c4_sh_nm.png', p_c4_sh, height = plotheight, width = plotwidth)
ggsave('plots/paper_figures/frequency_plots/cluster_frequencies/c5_sh_nm.pdf', p_c5_sh, height = plotheight, width = plotwidth)
ggsave('plots/frequency_plots/cluster_frequencies/c5_sh_nm.png', p_c5_sh, height = plotheight, width = plotwidth)

#######################
# save combined plots #
#######################

ggsave('plots/paper_figures/frequency_plots/cluster_frequencies/c1_c_nm_nolabels.pdf', p_c1_c, height = plotheight, width = plotwidth)
ggsave('plots/frequency_plots/cluster_frequencies/c1_c_nm_nolabels.png', p_c1_c, height = plotheight, width = plotwidth)
ggsave('plots/paper_figures/frequency_plots/cluster_frequencies/c2_c_nm_nolabels.pdf', p_c2_c, height = plotheight, width = plotwidth)
ggsave('plots/frequency_plots/cluster_frequencies/c2_c_nm_nolabels.png', p_c2_c, height = plotheight, width = plotwidth)
ggsave('plots/paper_figures/frequency_plots/cluster_frequencies/c3_c_nm_nolabels.pdf', p_c3_c, height = plotheight, width = plotwidth)
ggsave('plots/frequency_plots/cluster_frequencies/c3_c_nm_nolabels.png', p_c3_c, height = plotheight, width = plotwidth)
ggsave('plots/paper_figures/frequency_plots/cluster_frequencies/c4_c_nm_nolabels.pdf', p_c4_c, height = plotheight, width = plotwidth)
ggsave('plots/frequency_plots/cluster_frequencies/c4_c_nm_nolabels.png', p_c4_c, height = plotheight, width = plotwidth)
ggsave('plots/paper_figures/frequency_plots/cluster_frequencies/c5_c_nm_nolabels.pdf', p_c5_c, height = plotheight, width = plotwidth)
ggsave('plots/frequency_plots/cluster_frequencies/c5_c_nm_nolabels.png', p_c5_c, height = plotheight, width = plotwidth)

p_c1_c = p_c1_c + geom_text(aes(label=rownames(plotDF_c1),vjust=-0.4, size=2), show.legend = F)
p_c2_c = p_c2_c + geom_text(aes(label=rownames(plotDF_c2),vjust=-0.4, size=2), show.legend = F)
p_c3_c = p_c3_c + geom_text(aes(label=rownames(plotDF_c3),vjust=-0.4, size=2), show.legend = F)
p_c4_c = p_c4_c + geom_text(aes(label=rownames(plotDF_c4),vjust=-0.4, size=2), show.legend = F)
p_c5_c = p_c5_c + geom_text(aes(label=rownames(plotDF_c5),vjust=-0.4, size=2), show.legend = F)

ggsave('plots/paper_figures/frequency_plots/cluster_frequencies/c1_c_nm.pdf', p_c1_c, height = plotheight, width = plotwidth)
ggsave('plots/frequency_plots/cluster_frequencies/c1_c_nm.png', p_c1_c, height = plotheight, width = plotwidth)
ggsave('plots/paper_figures/frequency_plots/cluster_frequencies/c2_c_nm.pdf', p_c2_c, height = plotheight, width = plotwidth)
ggsave('plots/frequency_plots/cluster_frequencies/c2_c_nm.png', p_c2_c, height = plotheight, width = plotwidth)
ggsave('plots/paper_figures/frequency_plots/cluster_frequencies/c3_c_nm.pdf', p_c3_c, height = plotheight, width = plotwidth)
ggsave('plots/frequency_plots/cluster_frequencies/c3_c_nm.png', p_c3_c, height = plotheight, width = plotwidth)
ggsave('plots/paper_figures/frequency_plots/cluster_frequencies/c4_c_nm.pdf', p_c4_c, height = plotheight, width = plotwidth)
ggsave('plots/frequency_plots/cluster_frequencies/c4_c_nm.png', p_c4_c, height = plotheight, width = plotwidth)
ggsave('plots/paper_figures/frequency_plots/cluster_frequencies/c5_c_nm.pdf', p_c5_c, height = plotheight, width = plotwidth)
ggsave('plots/frequency_plots/cluster_frequencies/c5_c_nm.png', p_c5_c, height = plotheight, width = plotwidth)