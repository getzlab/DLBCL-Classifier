rm(list = ls())
source('src_R/load_libraries.R')

traininghistory = read.csv('model_training_history/combinedtraininghistory.txt', sep='\t')

meltedDF = melt(traininghistory, id.vars = c('PercentTrain'))

p<-ggplot(data=meltedDF, aes(x=PercentTrain, y=value, color=variable)) +
  geom_line() +
  xlab('Percent Trained') +
  ylab('Average Ensemble Loss') +
  scale_colour_discrete(name="Model Loss",
                      breaks=c("trainingloss", "validationloss"),
                      labels=c("Training", "Validation"))

pdf('plots/paper_figures/TrainingHistory.pdf')
print(p)
dev.off()

traininghistory$bin = cut(traininghistory$PercentTrain, breaks=seq(0, 1, 0.04))
traininghistory$PercentTrain_upper = as.numeric(traininghistory$bin) * 0.04
averageloss = aggregate(traininghistory$trainingloss, list(traininghistory$PercentTrain_upper), mean)
averagelossvalidation = aggregate(traininghistory$validationloss, list(traininghistory$PercentTrain_upper), mean)
colnames(averageloss) = c('PercentTrain', 'AverageTrainLoss')
averageloss$AverageValidationLoss = averagelossvalidation$x

meltedDF = melt(averageloss, id.vars = c('PercentTrain'))

p<-ggplot(data=meltedDF, aes(x=PercentTrain, y=value, color=variable)) +
  geom_line() + 
  geom_point() +
  xlab('Percent Trained') +
  ylab('Average Ensemble Loss (binned)') +
  scale_colour_discrete(name="Model Loss",
                        breaks=c("AverageTrainLoss", "AverageValidationLoss"),
                        labels=c("Training", "Validation"))

pdf('plots/paper_figures/TrainingHistory_binned.pdf')
print(p)
dev.off()

write.table(traininghistory, 'model_training_history/traininghistory_full.txt', sep='\t', quote=FALSE, row.names = FALSE)
write.table(averageloss, 'model_training_history/traininghistory_binned.txt', sep='\t', quote=FALSE, row.names=FALSE)

traininghistory_epochs = read.csv('model_training_history/combinedtrainhistory_epochs.txt', sep='\t')
std_devs = c(traininghistory_epochs$stdTrain, traininghistory_epochs$stdValidation)
std_errs = c(traininghistory_epochs$stdTrain / sqrt(traininghistory_epochs$num_models), 
             traininghistory_epochs$stdValidation / sqrt(traininghistory_epochs$num_models))
num_models = traininghistory_epochs$num_models

traininghistory_epochs = subset(traininghistory_epochs, select = -c(stdTrain, stdValidation, num_models))
meltedDF = melt(traininghistory_epochs, id.vars = c('Epoch'))
meltedDF$stdDev = std_devs
meltedDF$stdErr = std_errs

p<-ggplot(data=meltedDF, aes(x=Epoch, y=value, color=variable)) +
  geom_line() + 
  geom_point(size=0.5) +
  geom_errorbar(aes(ymin=value-stdErr, ymax=value+stdErr), width=0.0) +
  xlab('Epoch') +
  ylab('Average Ensemble Error') +
  scale_colour_discrete(name="Model Loss",
                        breaks=c("trainingloss", "validationloss"),
                        labels=c("Training", "Validation")) + 
  theme_bw() +
  theme(panel.grid = element_line(colour = "#C1C1C1")) +
  theme(axis.text.x = element_text(colour="grey20",size=10,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(size=10),  
        axis.title.x = element_text(colour="grey20",size=15,angle=0,hjust=.5,vjust=0,face="plain"),
        axis.title.y = element_text(size=15),
        legend.position = c(0.9, 0.9)) +
  scale_x_continuous(expand = c(0, 0), limits = c(-2.7, 45)) +
  scale_colour_manual(values = c("#497DE7","#E3140F"))

newDF = data.frame(traininghistory_epochs$Epoch, num_models)
colnames(newDF) = c('Epoch', 'Models')
p2 = ggplot(data=newDF, aes(x=Epoch, y=Models)) +
  geom_point() +
  theme_bw() + 
  theme(panel.grid = element_line(colour = "#C1C1C1")) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size=10),  
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=15),
        legend.position = c(0.9, 0.9),
        panel.grid.minor.y = element_blank()) +
  scale_x_continuous(expand = c(0, 0), limits = c(-2.7, 45))

p3 = plot_grid(p2,p, nrow=2, rel_widths = c(1,1), rel_heights = c(1, 4), align='v')


pdf('plots/paper_figures/TrainingHistory_epochs.pdf')
print(p)
dev.off()

jpeg('plots/paper_figures/TrainingHistory_epochs.jpeg')
print(p)
dev.off()

pdf('plots/paper_figures/TrainingHistory_epochs_modelsSurvived.pdf')
print(p3)
dev.off()
