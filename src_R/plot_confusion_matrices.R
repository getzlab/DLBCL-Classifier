rm(list = ls())
source("src_R/load_libraries.R")

preds_train = read.csv('evaluation_validation_set/confidence_adjusted_tables/NN_reducedV3.2_removeN5_nfeatures21_pMax0.94248563.tsv',
                       sep='\t', row.names=1)
new_labels = read.csv('data_tables/confidence_tables/baseline_probabilities.connectivity_based.sensitivity_power2.Aug_17_2022.tsv', 
                      sep='\t', row.names=1)

colnames(preds_train) = c('C1', 'C2', 'C3', 'C4', 'C5')
preds_train$Predicted_Cluster = apply(preds_train, 1, which.max)

Predicted_Label = paste('C', preds_train$Predicted_Cluster, sep='')
True_Label = paste('C', new_labels[rownames(preds_train), 'cluster'], sep='')
t_tab = table(Predicted_Label, True_Label)
plot_df = melt(t_tab)

p_train = ggplot(data =  plot_df, mapping = aes(x = Predicted_Label, y = True_Label)) +
          geom_tile(aes(fill = value), colour = "white") +
          geom_text(aes(label = sprintf("%1.0f", value)), vjust = 1) +
          scale_fill_gradient(low = "white", high = "green") +
          theme_bw() + theme(legend.position = "none") +
          scale_y_discrete(limits=rev) +
          ylab('True Cluster') + 
          xlab(paste('Predicted Cluster\n\nAccuracy = ', 
                     round(sum(Predicted_Label == True_Label) / length(Predicted_Label), 2) , 
                     sep=''))



ggsave('plots/paper_figures/confusion_matrix_train.pdf', p_train, height = 5, width = 5)
ggsave('plots/confusion_matrix_train.png', p_train, height = 5, width = 5)

preds_test = read.csv('evaluation_test_set/NN_reducedV3.2_removeN5_nfeatures21_testsetEval.tsv', sep='\t', row.names = 1)

Predicted_Label_Test = paste('C', preds_test$PredictedCluster, sep='')
True_Label_Test = paste('C', new_labels[rownames(preds_test), 'cluster'], sep='')
t_tab = table(Predicted_Label_Test, True_Label_Test)
plot_df = melt(t_tab)

p_test = ggplot(data =  plot_df, mapping = aes(x = Predicted_Label_Test, y = True_Label_Test)) +
  geom_tile(aes(fill = value), colour = "white") +
  geom_text(aes(label = sprintf("%1.0f", value)), vjust = 1) +
  scale_fill_gradient(low = "white", high = "green") +
  theme_bw() + theme(legend.position = "none") +
  scale_y_discrete(limits=rev) +
  ylab('True Cluster') + 
  xlab(paste('Predicted Cluster\n\nAccuracy = ', 
             round(sum(Predicted_Label_Test == True_Label_Test) / length(Predicted_Label_Test), 2) , 
             sep=''))

ggsave('plots/paper_figures/confusion_matrix_train_test.pdf', p_test, height = 5, width = 5)
ggsave('plots/confusion_matrix_train_test.png', p_test, height = 5, width = 5)
