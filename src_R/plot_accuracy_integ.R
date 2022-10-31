rm(list = ls())
source("src_R/load_libraries.R")

preds_train = read.csv('evaluation_validation_set/confidence_adjusted_tables/NN_reducedV3.4_removeN5_nfeatures21_pMax0.93856484.tsv',
                       sep='\t', row.names=1)
preds_test = read.csv('evaluation_test_set/NN_reducedV3.4_removeN5_nfeatures21_testsetEval.tsv', sep='\t', row.names=1)
labels = read.csv('data_tables/confidence_tables/baseline_probabilities.connectivity_based.sensitivity_power2.Sep_23_2022.tsv', 
                  sep='\t', row.names=1)

labels_train = labels[rownames(preds_train), ]

colnames(preds_train) = c('C1', 'C2', 'C3', 'C4', 'C5')

preds_train$Confidence = apply(preds_train, 1, max)
preds_train$PredictedCluster = apply(preds_train, 1, which.max)
preds_train$TrueCluster = labels_train$cluster

preds_train$Correctness = preds_train$TrueCluster == preds_train$PredictedCluster
preds_test$Correctness = preds_test$TrueCluster == preds_test$PredictedCluster

preds_all = rbind(preds_train, preds_test)

cutoffs = seq(0.20, 1, 0.01)

y_train = c()
y_test = c()
n_samples_train = c()
n_samples_test = c()

for(c in cutoffs){
  p_train = preds_train[preds_train$Confidence >= c, 'Correctness']
  num_c_train = sum(p_train)
  
  p_test = preds_test[preds_test$Confidence >= c, 'Correctness']
  num_c_test = sum(p_test)
  
  y_train = c(y_train, num_c_train / length(p_train))
  y_test = c(y_test, num_c_test / length(p_test))
  
  n_samples_train = c(n_samples_train, length(p_train))
  n_samples_test = c(n_samples_test, length(p_test))
}

plot_df_train = data.frame(cutoffs, y_train, n_samples_train)
colnames(plot_df_train) = c('ConfidenceCutoff', 'TotalAccuracy', 'Samples')

plot_df_test = data.frame(cutoffs, y_test, n_samples_test)
colnames(plot_df_test) = c('ConfidenceCutoff', 'TotalAccuracy', 'Samples')

p_train = ggplot(plot_df_train, aes(x=ConfidenceCutoff)) +
  geom_line(aes(y=TotalAccuracy * 550, color='red')) + 
  geom_line(aes(y=Samples)) +
  scale_y_continuous(name = 'Samples Above Threshold', sec.axis = sec_axis(~./550, name = 'Accuracy of Samples Above Threshold',
                                                                           breaks=seq(0,1,0.1)),
                     breaks=seq(0,550,55)) +
  theme_bw() +
  theme(axis.title.y.right = element_text(color = "red", size=18),
        axis.text.y.right = element_text(color = "red"),
        axis.ticks.y.right = element_line(color = "red"),
        axis.title.y.left = element_text(size=18),
        axis.title.x = element_text(size=18),
        legend.position = "none") +
  ggtitle('Train Set')


p_test = ggplot(plot_df_test, aes(x=ConfidenceCutoff)) +
  geom_line(aes(y=TotalAccuracy * 149, color='red')) + 
  geom_line(aes(y=Samples)) +
  scale_y_continuous(name = 'Samples Above Threshold', sec.axis = sec_axis(~./149, name = 'Accuracy of Samples Above Threshold')) +
  theme_bw() +
  theme(axis.title.y.right = element_text(color = "red", size=18),
         axis.text.y.right = element_text(color = "red"),
         axis.ticks.y.right = element_line(color = "red"),
         axis.title.y.left = element_text(size=18),
         axis.title.x = element_text(size=18),
         legend.position = "none") +
  ggtitle('Test Set')

ggsave('plots/threshold_conf_acc_train.png', p_train)
ggsave('plots/paper_figures/threshold_conf_acc_train.pdf', p_train)

ggsave('plots/test_set/threshold_conf_acc_test.png', p_test)
ggsave('plots/test_set/threshold_conf_acc_test.pdf', p_test)
