rm(list = ls())
source("src_R/load_libraries.R")
library(ggalluvial)
options(warn=-1)

new_labels = read.csv('data_tables/confidence_tables/baseline_probabilities.connectivity_based.sensitivity_power2.Sep_23_2022.tsv', 
                      sep='\t', row.names=1)
nm_labels = read.csv('data_tables/clustering_labels/NatMed.DLBCL.bestclus.txt', sep='\t', row.names=1)
preds_train = read.csv('evaluation_validation_set/confidence_adjusted_tables/NN_reducedV3.4_removeN5_nfeatures21_pMax0.93856484.tsv',
                       sep='\t', row.names=1)
preds_test = read.csv('evaluation_test_set/NN_reducedV3.4_removeN5_nfeatures21_testsetEval.tsv', sep='\t', row.names=1)
coo = read.csv('data_tables/phenotypes/COO_and_genetic_lables.txt', sep='\t', row.names=1)
cohorts = read.csv('data_tables/sample_sets/ShippStaudtSets.purity0.2.txt', row.names=1, sep='\t')
cohorts$dlbclass_cohort = 'Chapuy'
cohorts[cohorts$cohort != 'Shipp', 'dlbclass_cohort'] = 'Schmitz'

schmitz_s = rownames(cohorts)[cohorts$dlbclass_cohort == 'Schmitz']
chapuy_s = rownames(cohorts)[cohorts$dlbclass_cohort == 'Chapuy']

lymphgen = read.csv('data_tables/phenotypes/lymphgenclasses.tsv', sep='\t', row.names=1)
colnames(lymphgen)[4] = 'LymphGen_Class'

coo = coo[rownames(new_labels), ]
coo$COO = gsub('^Unclass$', 'UC', coo$COO)
coo$COO = gsub('^Unclassified$', 'UC', coo$COO)

colnames(preds_train) = c('C1','C2','C3','C4','C5')
preds_train$PredictedCluster = apply(preds_train,1,which.max)
preds_train$Confidence = apply(preds_train[, 1:5],1,max)
preds_train$TrueCluster = new_labels[rownames(preds_train), 'cluster']
preds_test$PredictedCluster = paste('C', preds_test$PredictedCluster, sep='') 

new_labels$cluster = paste('C', new_labels$cluster, sep='')
nm_labels$cluster = paste('C', nm_labels$cluster, sep='')

##############################
# NM clusters vs new clusters
##############################

nm_labels = nm_labels[rownames(nm_labels) %in% rownames(new_labels),]
new_labels_nmsub = new_labels[rownames(nm_labels), ]
tmp = data.frame(paste('C', preds_train$PredictedCluster, sep=''), preds_train$Confidence,
                 row.names = rownames(preds_train))
colnames(tmp) = c('PredictedCluster', 'Confidence')
all_preds = rbind(preds_test[, c('PredictedCluster', 'Confidence'), drop=F],
                  tmp)
all_preds = all_preds[rownames(nm_labels), ,drop=F]

nm_new_allu = data.frame(nm_labels$cluster, new_labels_nmsub$cluster, all_preds$PredictedCluster)
colnames(nm_new_allu) = c('NM_cluster', 'Combined_Clustering', 'DLBclass')

nm_new_allu$DLBclass_t = ifelse(all_preds$Confidence > 0.70, all_preds$PredictedCluster, 'Below\nThreshold')

nm_new_counts_2 = nm_new_allu %>% count(NM_cluster, Combined_Clustering, DLBclass, DLBclass_t)

nm_new_lodes_2 = to_lodes_form(nm_new_counts_2,
                             key = "ClusterSet",
                             axes = 1:(ncol(nm_new_counts_2) - 1))

tb = table(paste(nm_new_allu$Combined_Clustering, '_Combined  ', sep=''), paste(nm_new_allu$NM_cluster, '_NM', sep=''))
tb = melt(tb)

tb_2 = table(paste(nm_new_allu$DLBclass, '_DLBclass   ', sep=''), paste(nm_new_allu$Combined_Clustering, '_Combined', sep=''))
tb_2 = melt(tb_2)

tb_3 = table(paste(nm_new_allu$DLBclass, '_DLBclass   ', sep=''), paste(nm_new_allu$NM_cluster, '_NM', sep=''))
tb_3 = melt(tb_3)

tb_4 = table(ifelse(grepl('Below', nm_new_allu$DLBclass_t), nm_new_allu$DLBclass_t, paste(nm_new_allu$DLBclass_t, '_DLBclass(t)', sep='')), 
             paste(nm_new_allu$Combined_Clustering, '_Combined', sep=''))
tb_4 = melt(tb_4)

tb_5 = table(ifelse(grepl('Below', nm_new_allu$DLBclass_t), nm_new_allu$DLBclass_t, paste(nm_new_allu$DLBclass_t, '_DLBclass(t)', sep='')), 
             paste(nm_new_allu$NM_cluster, '_NM', sep=''))
tb_5 = melt(tb_5)

p_nm_new_2 = ggplot(data = nm_new_lodes_2,
                  aes(x = ClusterSet, stratum = stratum, alluvium = alluvium,
                      y = n, label = stratum)) +
  geom_alluvium(aes(fill = stratum)) +
  geom_stratum(color = "black", size = 0.3,
               fill=c(rep(rev(c("#803e98", "#00a4d1", "#f29123", "#4d872d", "#dd2929")),
                        (ncol(nm_new_counts_2) - 2)),
                      rev(c("#803e98", "#00a4d1", "#f29123", "#4d872d", "#dd2929", "#808080"))),
               alpha=1.0) +
  geom_text(stat = "stratum") +
  theme_bw() +
  scale_fill_manual(values = c("#803e98", "#00a4d1", "#f29123", "#4d872d", "#dd2929", "#808080"), name="Cluster Label") +
  theme(axis.text.y=element_blank(),
        axis.text.x = element_text(size=20),
        axis.ticks.y=element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5, size=20)
  ) +
  geom_text(stat = "stratum", aes(label = after_stat(count), vjust='top')) +
  ggtitle(paste('NM vs Combined_Clustering (n=', nrow(nm_labels) ,')', sep=''))

p_nm_new_2_table = ggplot(data = tb, mapping = aes(x = Var2, y = Var1)) +
  geom_tile(aes(fill = value)) +
  geom_text(aes(label = value), color="white" , vjust = 1) +
  theme_bw() + theme(legend.position = "none") + 
  scale_y_discrete(limits=rev) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank())

p_nm_new_2_table_2 = ggplot(data = tb_2, mapping = aes(x = Var2, y = Var1)) +
  geom_tile(aes(fill = value)) +
  geom_text(aes(label = value), color="white" , vjust = 1) +
  theme_bw() + theme(legend.position = "none") + 
  scale_y_discrete(limits=rev) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank())

p_nm_new_2_table_3 = ggplot(data = tb_3, mapping = aes(x = Var2, y = Var1)) +
  geom_tile(aes(fill = value)) +
  geom_text(aes(label = value), color="white" , vjust = 1) +
  theme_bw() + theme(legend.position = "none") + 
  scale_y_discrete(limits=rev) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank())

p_nm_new_2_table_4 = ggplot(data = tb_4, mapping = aes(x = Var2, y = Var1)) +
  geom_tile(aes(fill = value)) +
  geom_text(aes(label = value), color="white" , vjust = 1) +
  theme_bw() + theme(legend.position = "none") + 
  scale_y_discrete(limits=rev) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank())

p_nm_new_2_table_5 = ggplot(data = tb_5, mapping = aes(x = Var2, y = Var1)) +
  geom_tile(aes(fill = value)) +
  geom_text(aes(label = value), color="white" , vjust = 1) +
  theme_bw() + theme(legend.position = "none") + 
  scale_y_discrete(limits=rev) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank())

comb_2 = plot_grid(p_nm_new_2_table, p_nm_new_2_table_2,
                   p_nm_new_2_table_3, p_nm_new_2_table_4, p_nm_new_2_table_5,
                   rel_heights = c(0.20, 0.20, 0.20, 0.20, 0.20),
                   nrow = 5)

ggsave('plots/nm_newclus_alluvial_2.png', p_nm_new_2, width = 12, height=10)
ggsave('plots/paper_figures/nm_newclus_alluvial_2.pdf', p_nm_new_2, width = 12, height=10)

ggsave('plots/nm_newclus_confmats_2.png', comb_2, width = 6, height=10)
ggsave('plots/paper_figures/nm_newclus_confmats_2.pdf', comb_2, width = 6, height=10)

######################################
# clusters, init pred, init pred (t) #
######################################

tmp = data.frame(paste('C', preds_train$PredictedCluster, sep=''), preds_train$Confidence,
                 row.names = rownames(preds_train))
colnames(tmp) = c('PredictedCluster', 'Confidence')
all_preds = rbind(preds_test[, c('PredictedCluster', 'Confidence'), drop=F],
                  tmp)
all_preds = all_preds[rownames(new_labels), ,drop=F]

nm_new_allu = data.frame(new_labels$cluster, all_preds$PredictedCluster)
colnames(nm_new_allu) = c('Combined_Clustering', 'DLBclass')

nm_new_allu$DLBclass_t = ifelse(all_preds$Confidence > 0.70, all_preds$PredictedCluster, 'Below\nThreshold')

nm_new_counts_2 = nm_new_allu %>% count(Combined_Clustering, DLBclass, DLBclass_t)

nm_new_lodes_2 = to_lodes_form(nm_new_counts_2,
                               key = "ClusterSet",
                               axes = 1:(ncol(nm_new_counts_2) - 1))


tb_2 = table(paste(nm_new_allu$DLBclass, '_DLBclass   ', sep=''), paste(nm_new_allu$Combined_Clustering, '_Combined', sep=''))
tb_2 = melt(tb_2)

tb_4 = table(ifelse(grepl('Below', nm_new_allu$DLBclass_t), nm_new_allu$DLBclass_t, paste(nm_new_allu$DLBclass_t, '_DLBclass(t)', sep='')), 
             paste(nm_new_allu$Combined_Clustering, '_Combined', sep=''))
tb_4 = melt(tb_4)

p_nm_new_2 = ggplot(data = nm_new_lodes_2,
                    aes(x = ClusterSet, stratum = stratum, alluvium = alluvium,
                        y = n, label = stratum)) +
  geom_alluvium(aes(fill = stratum)) +
  geom_stratum(color = "black", size = 0.3,
               fill=c(rep(rev(c("#803e98", "#00a4d1", "#f29123", "#4d872d", "#dd2929")),
                          (ncol(nm_new_counts_2) - 2)),
                      rev(c("#803e98", "#00a4d1", "#f29123", "#4d872d", "#dd2929", "#808080"))),
               alpha=1.0) +
  geom_text(stat = "stratum") +
  theme_bw() +
  scale_fill_manual(values = c("#803e98", "#00a4d1", "#f29123", "#4d872d", "#dd2929", "#808080"), name="Cluster Label") +
  theme(axis.text.y=element_blank(),
        axis.text.x = element_text(size=20),
        axis.ticks.y=element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5, size=20)
  ) +
  geom_text(stat = "stratum", aes(label = after_stat(count), vjust='top')) +
  ggtitle(paste('Combined_Clustering vs DLBclass (n=', nrow(all_preds) ,')', sep=''))

p_nm_new_2_table_2 = ggplot(data = tb_2, mapping = aes(x = Var2, y = Var1)) +
  geom_tile(aes(fill = value)) +
  geom_text(aes(label = value), color="white" , vjust = 1) +
  theme_bw() + theme(legend.position = "none") + 
  scale_y_discrete(limits=rev) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank())

p_nm_new_2_table_4 = ggplot(data = tb_4, mapping = aes(x = Var2, y = Var1)) +
  geom_tile(aes(fill = value)) +
  geom_text(aes(label = value), color="white" , vjust = 1) +
  theme_bw() + theme(legend.position = "none") + 
  scale_y_discrete(limits=rev) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank())

comb_2 = plot_grid(p_nm_new_2_table_2, p_nm_new_2_table_4,
                   rel_heights = c(0.20, 0.20, 0.20, 0.20, 0.20),
                   nrow = 2)

ggsave('plots/newclus_dlbclass(t)_alluvial_2.png', p_nm_new_2, width = 12, height=10)
ggsave('plots/paper_figures/newclus_dlbclass(t)_alluvial_2.pdf', p_nm_new_2, width = 12, height=10)

ggsave('plots/newclus_dlbclass(t)_confmats_2.png', comb_2, width = 6, height=10)
ggsave('plots/paper_figures/newclus_dlbclass(t)_confmats_2.pdf', comb_2, width = 6, height=10)

##############################
# New clusters vs initial pred
##############################

cf_preds_test = preds_test
cf_preds_train = preds_train

cf_preds_test$PredictedCluster = paste(as.character(cf_preds_test$PredictedCluster),'_Test', sep='')
cf_preds_train$PredictedCluster = paste('C', as.character(cf_preds_train$PredictedCluster),'_Train', sep='')

#test_train_df = rbind(cf_preds_test[, c('PredictedCluster', 'TrueCluster')], cf_preds_train[, c('PredictedCluster', 'TrueCluster')])
train_df = cf_preds_train[, c('PredictedCluster', 'TrueCluster')]
test_df = cf_preds_test[, c('PredictedCluster', 'TrueCluster')]

colnames(train_df) = c('DLBclass', 'ClusterLabel')
colnames(test_df) = c('DLBclass', 'ClusterLabel')

train_df$ClusterLabel = paste('C', train_df$ClusterLabel, sep='')
test_df$ClusterLabel = paste('C', test_df$ClusterLabel, sep='')

train_counts = train_df %>% count(DLBclass, ClusterLabel)
test_counts = test_df %>% count(DLBclass, ClusterLabel)

train_lodes = to_lodes_form(train_counts,
                            key = "Xaxis",
                            axes = 1:(ncol(train_counts) - 1))

test_lodes = to_lodes_form(test_counts,
                           key = "Xaxis",
                           axes = 1:(ncol(test_counts) - 1))

tb = table(paste(train_df$DLBclass, '_DLBclass', sep=''), paste(train_df$ClusterLabel, '_ClusLabel', sep=''))
tb = melt(tb)

tb_test = table(paste(test_df$DLBclass, '_DLBclass', sep=''), paste(test_df$ClusterLabel, '_ClusLabel', sep=''))
tb_test = melt(tb_test)

p_train = ggplot(data = train_lodes,
                 aes(x = Xaxis, stratum = stratum, alluvium = alluvium,
                     y = n, label = stratum)) +
  geom_flow(aes(fill = stratum)) +
  scale_alpha_discrete(range=c(0.9, 0.4)) +
  geom_stratum(color = "black", size = 0.3,
               fill=c(rev(c("#803e98", "#00a4d1", "#f29123", "#4d872d", "#dd2929")),
                      rev(c("#803e98", "#00a4d1", "#f29123", "#4d872d", "#dd2929"))),
               alpha=1.0) +
  geom_text(stat = "stratum") +
  theme_bw() +
  scale_fill_manual(values = c("#803e98", "#00a4d1", "#f29123", "#4d872d", "#dd2929",
                               "#803e98", "#00a4d1", "#f29123", "#4d872d", "#dd2929"), name="Cluster Label") +
  theme(axis.text.y=element_blank(),
        axis.text.x = element_text(size=20),
        axis.ticks.y=element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5, size=20)
  ) +
  geom_text(stat = "stratum", aes(label = after_stat(count), vjust='top')) +
  ggtitle(paste('Cluster Labels vs DLBclass (Train Set n=', nrow(train_df) ,')', sep=''))

p_test = ggplot(data = test_lodes,
                 aes(x = Xaxis, stratum = stratum, alluvium = alluvium,
                     y = n, label = stratum)) +
  geom_flow(aes(fill = stratum)) +
  scale_alpha_discrete(range=c(0.9, 0.4)) +
  geom_stratum(color = "black", size = 0.3,
               fill=c(rev(c("#803e98", "#00a4d1", "#f29123", "#4d872d", "#dd2929")),
                      rev(c("#803e98", "#00a4d1", "#f29123", "#4d872d", "#dd2929"))),
               alpha=1.0) +
  geom_text(stat = "stratum") +
  theme_bw() +
  scale_fill_manual(values = c("#803e98", "#00a4d1", "#f29123", "#4d872d", "#dd2929",
                               "#803e98", "#00a4d1", "#f29123", "#4d872d", "#dd2929"), name="Cluster Label") +
  theme(axis.text.y=element_blank(),
        axis.text.x = element_text(size=20),
        axis.ticks.y=element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5, size=20)
  ) +
  geom_text(stat = "stratum", aes(label = after_stat(count), vjust='top')) +
  ggtitle(paste('Cluster Labels vs DLBclass (Test Set n=', nrow(test_df) ,')', sep=''))

p_train_table = ggplot(data = tb, mapping = aes(x = Var1, y = Var2)) +
  geom_tile(aes(fill = value)) +
  geom_text(aes(label = value), color="white" , vjust = 1) +
  theme_bw() + theme(legend.position = "none",
                     axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  scale_y_discrete(limits=rev)

p_test_table = ggplot(data = tb_test, mapping = aes(x = Var1, y = Var2)) +
  geom_tile(aes(fill = value)) +
  geom_text(aes(label = value), color="white" , vjust = 1) +
  theme_bw() + theme(legend.position = "none",
                     axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  scale_y_discrete(limits=rev)

comb_train = plot_grid(p_train, p_train_table, rel_heights = c(0.50, 0.50), nrow = 2)
comb_test = plot_grid(p_test, p_test_table, rel_heights = c(0.50, 0.50), nrow = 2)

ggsave('plots/newclus_initpred_train_alluvial.png', comb_train, width = 12, height=8)
ggsave('plots/paper_figures/newclus_initpred_train_alluvial.pdf', comb_train, width = 12, height=8)

ggsave('plots/newclus_initpred_test_alluvial.png', comb_test, width = 12, height=8)
ggsave('plots/paper_figures/newclus_initpred_test_alluvial.pdf', comb_test, width = 12, height=8)

##############################
# initial pred vs coo
##############################

cf_preds_test = preds_test
cf_preds_train = preds_train

cf_preds_test$coo = coo[rownames(cf_preds_test), 'COO']
cf_preds_train$coo = coo[rownames(cf_preds_train), 'COO']

cf_preds_test = cf_preds_test[cf_preds_test$coo != "na", ]
cf_preds_train = cf_preds_train[cf_preds_train$coo != "na", ]

cf_preds_test$HC = ifelse(cf_preds_test$Confidence > 0.70, 'HC', 'LC')
cf_preds_train$HC = ifelse(cf_preds_train$Confidence > 0.70, 'HC', 'LC')

cf_preds_test$TrueCluster = paste('C', cf_preds_test$TrueCluster, sep='')
cf_preds_train$TrueCluster = paste('C', cf_preds_train$TrueCluster, sep='')

cf_preds_train$PredictedCluster = paste('C', cf_preds_train$PredictedCluster, sep='')
cf_preds_test$PredictedCluster = paste(cf_preds_test$PredictedCluster, sep='')

cf_preds_test_conf = cf_preds_test
cf_preds_train_conf = cf_preds_train

cf_preds_train_conf$PredictedCluster = paste(cf_preds_train_conf$PredictedCluster, '_', cf_preds_train_conf$HC, sep='')
cf_preds_test_conf$PredictedCluster = paste(cf_preds_test_conf$PredictedCluster, '_', cf_preds_test_conf$HC, sep='')

conf_coo_df = rbind(cf_preds_train_conf[, c('PredictedCluster', 'coo', 'HC')], cf_preds_test_conf[, c('PredictedCluster', 'coo', 'HC')])

test_train_df = rbind(cf_preds_test[, c('PredictedCluster', 'coo')], cf_preds_train[, c('PredictedCluster', 'coo')])
test_train_df_2 = test_train_df
test_train_df_2$alpha = 'Schmitz'
test_train_df_2[rownames(test_train_df_2) %in% chapuy_s, 'alpha'] = 'Chapuy'
test_train_df_2$PredictedCluster = paste(test_train_df_2$PredictedCluster, '_', test_train_df_2$alpha, sep='')
colnames(test_train_df) = c('DLBclass', 'COO')
colnames(test_train_df_2) = c('DLBclass', 'COO', 'alpha')
colnames(conf_coo_df) = c('DLBclass', 'COO', 'alpha')

#test_train_df$Set = gsub('C\\w_', '', test_train_df$ClusterLabel)
#test_train_counts = test_train_df %>% count(ClusterLabel, DLBclass, Set)

test_train_counts = test_train_df %>% count(DLBclass, COO)
test_train_counts_2 = test_train_df %>% count(COO, DLBclass)
test_train_counts_3 = test_train_df_2 %>% count(DLBclass, COO)
conf_counts = conf_coo_df %>% count(DLBclass, COO)

test_train_lodes = to_lodes_form(test_train_counts,
                                 key = "Xaxis",
                                 axes = 1:(ncol(test_train_counts) - 1))
test_train_lodes_2 = to_lodes_form(test_train_counts_2,
                                   key = "Xaxis",
                                   axes = 1:(ncol(test_train_counts_2) - 1))
test_train_lodes_3 = to_lodes_form(test_train_counts_3,
                                   key = "Xaxis",
                                   axes = 1:(ncol(test_train_counts_3) - 1))
conf_lodes = to_lodes_form(conf_counts,
                           key = "Xaxis",
                           axes = 1:(ncol(conf_counts) - 1))

test_train_lodes$alpha = gsub('C\\w_', '', test_train_lodes$stratum)
test_train_lodes[30:nrow(test_train_lodes), 'alpha'] = 'Train'

#test_train_lodes_2$alpha = gsub('C\\w_', '', test_train_lodes_2$stratum)
test_train_lodes_2$alpha = "Train"

test_train_lodes_3$alpha = gsub('C\\w_', '', test_train_lodes_3$stratum)
test_train_lodes_3[test_train_lodes_3$Xaxis == 'COO', 'alpha'] = 'none'

conf_lodes$alpha = gsub('C\\w_', '', conf_lodes$stratum)
conf_lodes[conf_lodes$Xaxis == 'COO', 'alpha'] = 'none'

tb = table(paste(test_train_df$DLBclass, '_DLBclass', sep=''), test_train_df$COO)
tb = melt(tb)

p_coo = ggplot(data = test_train_lodes,
               aes(x = Xaxis, stratum = stratum, alluvium = alluvium,
                   y = n, label = stratum, alpha=alpha)) +
  geom_flow(aes(fill = stratum)) +
  scale_alpha_discrete(range=c(0.9, 0.4)) +
  geom_stratum(color = "black", size = 0.3,
               fill=c(rev(c("#803e98", "#00a4d1", "#f29123", "#4d872d", "#dd2929")),
                      rev(c("#fca59d", "#194bff" ,"#f5fc72"))
                      ),
               alpha=1.0) +
  geom_text(stat = "stratum") +
  theme_bw() +
  scale_fill_manual(values = c("#803e98", "#00a4d1", "#f29123", "#4d872d", "#dd2929", 
                               "#000000","#000000","#000000"), name="Cluster Label") +
  theme(axis.text.y=element_blank(),
        axis.text.x = element_text(size=20),
        axis.ticks.y=element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5, size=20)
  ) +
  geom_text(stat = "stratum", aes(label = after_stat(count), vjust='top')) +
  ggtitle(paste('DLBclass vs COO (n=', nrow(test_train_df),')', sep=''))

p_coo_2 = ggplot(data = test_train_lodes_2,
               aes(x = Xaxis, stratum = stratum, alluvium = alluvium,
                   y = n, label = stratum, alpha=alpha)) +
  scale_alpha_discrete(range=c(0.8, 1.0)) +
  geom_flow(aes(fill = stratum)) +
  geom_stratum(color = "black", size = 0.3,
               fill=c(rev(c("#fca59d", "#194bff", "#f5fc72")),
                      rev(c("#803e98", "#00a4d1", "#f29123", "#4d872d", "#dd2929"))
                      ),
               alpha=1.0) +
  geom_text(stat = "stratum") +
  theme_bw() +
  scale_fill_manual(values = c("#fca59d", "#194bff", "#f5fc72",
                               "#803e98", "#00a4d1", "#f29123", "#4d872d", "#dd2929"), name="Cluster Label") +
  theme(axis.text.y=element_blank(),
        axis.text.x = element_text(size=20),
        axis.ticks.y=element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5, size=20)
  ) +
  geom_text(stat = "stratum", aes(label = after_stat(count), vjust='top')) +
  ggtitle(paste('DLBclass vs COO (n=', nrow(test_train_df) ,')', sep=''))

p_coo_2_table = ggplot(data = tb, mapping = aes(x = Var1, y = Var2)) +
  geom_tile(aes(fill = value)) +
  geom_text(aes(label = value), color="white" , vjust = 1) +
  theme_bw() + theme(legend.position = "none",
                     axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  scale_y_discrete(limits=rev)

p_coo_3 = ggplot(data = test_train_lodes_3,
              aes(x = Xaxis, stratum = stratum, alluvium = alluvium,
                  y = n, label = stratum, alpha=alpha)) +
  geom_flow(aes(fill = stratum)) +
  scale_alpha_discrete(range=c(0.9, 0.3, 1.0)) +
  geom_stratum(color = "black", size = 0.3,
               fill = c(rev(c("#803e98", "#803e98",
                              "#00a4d1", "#00a4d1",
                              "#f29123", "#f29123", 
                              "#4d872d", "#4d872d",
                              "#dd2929", "#dd2929")), 
                        rev(c("#fca59d", "#194bff", "#f5fc72"))),
               alpha=c(rep(c(0.5, 1.0), 5), rep(1.0, 3))) +
  geom_text(stat = "stratum") +
  scale_fill_manual(values = c(c("#803e98", "#803e98",
                                 "#00a4d1", "#00a4d1",
                                 "#f29123", "#f29123",
                                 "#4d872d", "#4d872d",
                                 "#dd2929", "#dd2929"), 
                               rep("#FFFFFF", 3)), name="Cluster Label") +
  theme_bw() +
  theme(axis.text.y=element_blank(),
        axis.text.x = element_text(size=20),
        axis.ticks.y=element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5, size=20)
  ) +
  geom_text(stat = "stratum", aes(label = after_stat(count), vjust='top')) +
  ggtitle(paste('DLBclass (Cohort) vs COO Predictions (n=', nrow(test_train_df_2) ,')', sep=''))

p_coo_conf = ggplot(data = conf_lodes,
                 aes(x = Xaxis, stratum = stratum, alluvium = alluvium,
                     y = n, label = stratum, alpha=alpha)) +
  geom_flow(aes(fill = stratum)) +
  scale_alpha_discrete(range=c(0.9, 0.3, 1.0)) +
  geom_stratum(color = "black", size = 0.3,
               fill = c(rev(c("#803e98", "#803e98",
                              "#00a4d1", "#00a4d1",
                              "#f29123", "#f29123", 
                              "#4d872d", "#4d872d",
                              "#dd2929", "#dd2929")), 
                        rev(c("#fca59d", "#194bff", "#f5fc72"))),
               alpha=c(rep(c(0.5, 1.0), 5), rep(1.0, 3))) +
  geom_text(stat = "stratum") +
  scale_fill_manual(values = c(c("#803e98", "#803e98",
                                 "#00a4d1", "#00a4d1",
                                 "#f29123", "#f29123",
                                 "#4d872d", "#4d872d",
                                 "#dd2929", "#dd2929"), 
                               rep("#FFFFFF", 3)), name="Cluster Label") +
  theme_bw() +
  theme(axis.text.y=element_blank(),
        axis.text.x = element_text(size=20),
        axis.ticks.y=element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5, size=20)
  ) +
  geom_text(stat = "stratum", aes(label = after_stat(count), vjust='top')) +
  ggtitle(paste('DLBclass (Confidence) vs COO Predictions (n=', nrow(conf_coo_df) ,')', sep=''))

comb_2 = plot_grid(p_coo_2, p_coo_2_table, rel_heights = c(0.50, 0.50), nrow = 2)

ggsave('plots/initpred_coo_alluvial.png', p_coo, width = 12, height=8)
ggsave('plots/paper_figures/initpred_coo_alluvial.pdf', p_coo, width = 12, height=8)

ggsave('plots/initpred_coo_alluvial_2.png', comb_2, width = 12, height=8)
ggsave('plots/paper_figures/initpred_coo_alluvial_2.pdf', comb_2, width = 12, height=8)

ggsave('plots/initpred_coo_alluvial_cohort.png', p_coo_3, width = 12, height=8)
ggsave('plots/paper_figures/initpred_coo_alluvial_cohort.pdf', p_coo_3, width = 12, height=8)

ggsave('plots/initpred_coo_alluvial_conf.png', p_coo_conf, width = 12, height=8)
ggsave('plots/paper_figures/initpred_coo_alluvial_conf.pdf', p_coo_conf, width = 12, height=8)

######################################
# initial pred vs true cluster w/ conf
######################################

cf_preds_test = preds_test
cf_preds_train = preds_train

cf_preds_test$HC = ifelse(cf_preds_test$Confidence > 0.80, 'HC', 'LC')
cf_preds_train$HC = ifelse(cf_preds_train$Confidence > 0.80, 'HC', 'LC')

cf_preds_train$PredictedCluster = paste('C', cf_preds_train$PredictedCluster, '_', cf_preds_train$HC, sep='')
cf_preds_test$PredictedCluster = paste(cf_preds_test$PredictedCluster, '_', cf_preds_test$HC, sep='')

cf_preds_test$TrueCluster = paste('C', cf_preds_test$TrueCluster, sep='')
cf_preds_train$TrueCluster = paste('C', cf_preds_train$TrueCluster, sep='')

test_train_df = rbind(cf_preds_test[, c('PredictedCluster', 'TrueCluster')], cf_preds_train[, c('PredictedCluster', 'TrueCluster')])
colnames(test_train_df) = c('DLBclass', 'ClusterLabel')

#test_train_df$Set = gsub('C\\w_', '', test_train_df$ClusterLabel)
#test_train_counts = test_train_df %>% count(ClusterLabel, DLBclass, Set)

test_train_counts = test_train_df %>% count(DLBclass, ClusterLabel)
test_train_counts_2 = test_train_df %>% count(ClusterLabel, DLBclass)

test_train_lodes = to_lodes_form(test_train_counts,
                                 key = "Xaxis",
                                 axes = 1:(ncol(test_train_counts) - 1))

test_train_lodes_2 = to_lodes_form(test_train_counts_2,
                                 key = "Xaxis",
                                 axes = 1:(ncol(test_train_counts_2) - 1))

test_train_lodes$alpha = gsub('C\\w_', '', test_train_lodes$stratum)
test_train_lodes[38:nrow(test_train_lodes), 'alpha'] = 'Train'

tb = table(paste(test_train_df$DLBclass, '_DLBclass', sep=''), paste(test_train_df$ClusterLabel, '_ClusLabel', sep=''))
tb = melt(tb)

p_hc_lc = ggplot(data = test_train_lodes,
                 aes(x = Xaxis, stratum = stratum, alluvium = alluvium,
                     y = n, label = stratum, alpha=alpha)) +
  geom_flow(aes(fill = stratum)) +
  scale_alpha_discrete(range=c(0.9, 0.3)) +
  geom_stratum(color = "black", size = 0.3,
               fill=c(rev(c("#803e98", "#803e98", "#00a4d1", "#00a4d1", "#f29123", "#f29123", "#4d872d", "#4d872d", "#dd2929", "#dd2929")),
                      rev(c("#803e98", "#00a4d1", "#f29123", "#4d872d", "#dd2929"))),
               alpha=1.0) +
  geom_text(stat = "stratum") +
  theme_bw() +
  scale_fill_manual(values = c("#803e98", "#803e98", "#00a4d1", "#00a4d1", "#f29123", "#f29123", 
                               "#4d872d", "#4d872d", "#dd2929", "#dd2929",
                               "#803e98", "#00a4d1", "#f29123", "#4d872d", "#dd2929"), name="Cluster Label") +
  theme(axis.text.y=element_blank(),
        axis.text.x = element_text(size=20),
        axis.ticks.y=element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5, size=20)
  ) +
  geom_text(stat = "stratum", aes(label = after_stat(count), vjust='top')) +
  ggtitle('Cluster Labels vs DLBclass')

p_hc_lc_2 = ggplot(data = test_train_lodes_2,
                 aes(x = Xaxis, stratum = stratum, alluvium = alluvium,
                     y = n, label = stratum)) +
  geom_flow(aes(fill = stratum)) +
  geom_stratum(color = "black", size = 0.3,
               fill=c(rev(c("#803e98", "#00a4d1", "#f29123", "#4d872d", "#dd2929")),
                        rev(c("#803e98", "#803e98", "#00a4d1", "#00a4d1", "#f29123", "#f29123", "#4d872d", "#4d872d", "#dd2929", "#dd2929"))),
               alpha=1.0) +
  geom_text(stat = "stratum") +
  theme_bw() +
  scale_fill_manual(values = c("#803e98", "#00a4d1", "#f29123", "#4d872d", "#dd2929",
                               "#803e98", "#803e98", "#00a4d1", "#00a4d1", "#f29123", "#f29123", "#4d872d", "#4d872d", "#dd2929", "#dd2929"),
                    name="Cluster Label") +
  theme(axis.text.y=element_blank(),
        axis.text.x = element_text(size=20),
        axis.ticks.y=element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5, size=20)
  ) +
  geom_text(stat = "stratum", aes(label = after_stat(count), vjust='top')) +
  ggtitle(paste('Cluster Labels vs DLBclass (n=', nrow(test_train_df) ,')', sep=''))

p_hc_lc_2_table = ggplot(data = tb, mapping = aes(x = Var1, y = Var2)) +
  geom_tile(aes(fill = value)) +
  geom_text(aes(label = value), color="white" , vjust = 1) +
  theme_bw() + theme(legend.position = "none",
                     axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  scale_y_discrete(limits=rev)

comb_2 = plot_grid(p_hc_lc_2, p_hc_lc_2_table, rel_heights = c(0.50, 0.50), nrow = 2)

ggsave('plots/initpred_newclus_hclc_alluvial.png', p_hc_lc, width = 12, height=8)
ggsave('plots/paper_figures/initpred_newclus_hclc_alluvial.pdf', p_hc_lc, width = 12, height=8)

ggsave('plots/initpred_newclus_hclc_alluvial_2.png', comb_2, width = 12, height=8)
ggsave('plots/paper_figures/initpred_newclus_hclc_alluvial_2.pdf', comb_2, width = 12, height=8)

######################################
# initial pred vs lymphgen
######################################

HC_THRESH = 0.70

lg_preds_train = preds_train
lg_preds_train$HC = ifelse(lg_preds_train$Confidence > HC_THRESH, 'HC', 'LC')
lg_preds_train$PredictedCluster = paste('C', lg_preds_train$PredictedCluster, sep='')

lg_preds_test = preds_test
lg_preds_test$HC = ifelse(lg_preds_test$Confidence > HC_THRESH, 'HC', 'LC')
#lg_preds_test$PredictedCluster = paste(lg_preds_test$PredictedCluster, '_', lg_preds_test$HC, sep='')

all_preds = rbind(lg_preds_test[,c('PredictedCluster', 'HC'),drop=FALSE], lg_preds_train[, c('PredictedCluster', 'HC'), drop=FALSE])
lymphgen = lymphgen[rownames(all_preds), ]

preds_lg_df = data.frame(all_preds$PredictedCluster)
preds_lg_df$LymphgenClass = lymphgen$LymphGen_Class
preds_lg_df$HC = all_preds$HC
colnames(preds_lg_df) = c('DLBclass', 'Lymphgen_Class', 'HC')
rownames(preds_lg_df) = rownames(all_preds)
preds_lg_df[grepl('/', preds_lg_df$Lymphgen_Class), 'Lymphgen_Class'] = 'Ambiguous' 
preds_lg_df$Lymphgen_Class = factor(preds_lg_df$Lymphgen_Class, levels = c('BN2', 'A53', 'EZB', 'ST2', 'MCD', 'N1', 'Ambiguous', 'Other'))
preds_lg_df = preds_lg_df[order(preds_lg_df$Lymphgen_Class), ]

preds_lg_df_hc = preds_lg_df[preds_lg_df$HC == 'HC', ]
preds_lg_df_lc = preds_lg_df[preds_lg_df$HC == 'LC', ]

preds_lg_df$DLBclass = paste(preds_lg_df$DLBclass, '_', preds_lg_df$HC, sep='')

preds_lg_counts = preds_lg_df %>% count(DLBclass, Lymphgen_Class)
preds_lg_counts_2 = preds_lg_df_hc %>% count(DLBclass, Lymphgen_Class)
preds_lg_counts_3 = preds_lg_df_lc %>% count(DLBclass, Lymphgen_Class)

preds_lg_lodes = to_lodes_form(preds_lg_counts,
                                 key = "Xaxis",
                                 axes = 1:(ncol(preds_lg_counts) - 1))

preds_lg_lodes_2 = to_lodes_form(preds_lg_counts_2,
                                 key = "Xaxis",
                                 axes = 1:(ncol(preds_lg_counts_2) - 1))

preds_lg_lodes_3 = to_lodes_form(preds_lg_counts_3,
                                 key = "Xaxis",
                                 axes = 1:(ncol(preds_lg_counts_3) - 1))

tb = table(paste(preds_lg_df_hc$DLBclass, '_DLBclass', sep=''), preds_lg_df_hc$Lymphgen_Class)
tb = melt(tb)

preds_lg_lodes$alpha = gsub('C\\w_', '', preds_lg_lodes$stratum)
preds_lg_lodes[preds_lg_lodes$Xaxis == 'Lymphgen_Class', 'alpha'] = 'none'

p_lg = ggplot(data = preds_lg_lodes,
                 aes(x = Xaxis, stratum = stratum, alluvium = alluvium,
                     y = n, label = stratum, alpha=alpha)) +
  geom_flow(aes(fill = stratum)) +
  scale_alpha_discrete(range=c(0.9, 0.3, 1.0)) +
  geom_stratum(color = "black", size = 0.3,
               fill = c(rev(c("#803e98", "#803e98",
                              "#00a4d1", "#00a4d1",
                              "#f29123", "#f29123", 
                              "#4d872d", "#4d872d",
                              "#dd2929", "#dd2929")), 
                        rep("#FFFFFF", 8)),
               alpha=c(rep(c(0.5, 1.0), 5), rep(1.0, 8))) +
  geom_text(stat = "stratum") +
  scale_fill_manual(values = c(c("#803e98", "#803e98",
                                 "#00a4d1", "#00a4d1",
                                 "#f29123", "#f29123",
                                 "#4d872d", "#4d872d",
                                 "#dd2929", "#dd2929"), 
                               rep("#FFFFFF", 20)), name="Cluster Label") +
  theme_bw() +
  theme(axis.text.y=element_blank(),
        axis.text.x = element_text(size=20),
        axis.ticks.y=element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5, size=20)
  ) +
  geom_text(stat = "stratum", aes(label = after_stat(count), vjust='top')) +
  ggtitle(paste('DLBclass vs LymphGen Predictions (n=', nrow(preds_lg_df) ,')', sep=''))

p_lg_hc_only = ggplot(data = preds_lg_lodes_2,
                      aes(x = Xaxis, stratum = stratum, alluvium = alluvium,
                          y = n, label = stratum)) +
  geom_flow(aes(fill = stratum)) +
  geom_stratum(color = "black", size = 0.3,
               fill = c(rev(c("#803e98", 
                              "#00a4d1", 
                              "#f29123", 
                              "#4d872d",
                              "#dd2929")), 
                        rep("#FFFFFF", 8)),
               alpha=1.0) +
  geom_text(stat = "stratum") +
  scale_fill_manual(values = c(c("#803e98", 
                                 "#00a4d1", 
                                 "#f29123", 
                                 "#4d872d",
                                 "#dd2929"), 
                               rep("#FFFFFF", 8)), name="Cluster Label") +
  theme_bw() +
  theme(axis.text.y=element_blank(),
        axis.text.x = element_text(size=20),
        axis.ticks.y=element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5, size=20)
  ) +
  geom_text(stat = "stratum", aes(label = after_stat(count), vjust='top')) +
  ggtitle(paste('DLBclass vs LymphGen Predictions (HC>=0.7 n=', nrow(preds_lg_df_hc) ,')', sep=''))

p_lg_lc_only = ggplot(data = preds_lg_lodes_3,
                      aes(x = Xaxis, stratum = stratum, alluvium = alluvium,
                          y = n, label = stratum)) +
  geom_flow(aes(fill = stratum)) +
  geom_stratum(color = "black", size = 0.3,
               fill = c(rev(c("#803e98", 
                              "#00a4d1", 
                              "#f29123", 
                              "#4d872d",
                              "#dd2929")), 
                        rep("#FFFFFF", 8)),
               alpha=1.0) +
  geom_text(stat = "stratum") +
  scale_fill_manual(values = c(c("#803e98", 
                                 "#00a4d1", 
                                 "#f29123", 
                                 "#4d872d",
                                 "#dd2929"), 
                               rep("#FFFFFF", 8)), name="Cluster Label") +
  theme_bw() +
  theme(axis.text.y=element_blank(),
        axis.text.x = element_text(size=20),
        axis.ticks.y=element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5, size=20)
  ) +
  geom_text(stat = "stratum", aes(label = after_stat(count), vjust='top')) +
  ggtitle(paste('DLBclass vs LymphGen Predictions (LC<0.7 n=', nrow(preds_lg_df_lc) ,')', sep=''))

p_lg_hc_only_table = ggplot(data = tb, mapping = aes(x = Var1, y = Var2)) +
  geom_tile(aes(fill = value)) +
  geom_text(aes(label = value), color="white" , vjust = 1) +
  theme_bw() + theme(legend.position = "none",
                     axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  scale_y_discrete(limits=rev)

comb_2 = plot_grid(p_lg_hc_only, p_lg_hc_only_table, rel_heights = c(0.50, 0.50), nrow = 2)

ggsave(paste('plots/initpred_lg_alluvial_hconly_', HC_THRESH ,'.png', sep=''), comb_2, width =10, height=10)
ggsave(paste('plots/paper_figures/initpred_lg_alluvial_hconly_', HC_THRESH ,'.pdf', sep=''), comb_2, width = 10, height=10)

ggsave(paste('plots/initpred_lg_alluvial_hconly_', HC_THRESH ,'_notable.png', sep=''), p_lg_hc_only, width = 10, height=10)
ggsave(paste('plots/paper_figures/initpred_lg_alluvial_hconly_', HC_THRESH ,'_notable.pdf', sep=''), p_lg_hc_only, width = 10, height=10)

ggsave(paste('plots/initpred_lg_alluvial_lconly_', HC_THRESH ,'_notable.png', sep=''), p_lg_lc_only, width = 10, height=10)
ggsave(paste('plots/paper_figures/initpred_lg_alluvial_lconly_', HC_THRESH ,'_notable.pdf', sep=''), p_lg_lc_only, width = 10, height=10)

ggsave('plots/initpred_lg_alluvial.png', p_lg, width = 10, height=10)
ggsave('plots/paper_figures/initpred_lg_alluvial.pdf', p_lg, width = 10, height=10)

