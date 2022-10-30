rm(list = ls())
source('src_R/load_libraries.R')

performancetable = read.csv('evaluation_validation_set/allmodelsevaluated.tsv', sep='\t')
rownames(performancetable) = performancetable$experiment

performancetable$experiment <- 
  factor(performancetable$experiment, levels = performancetable$experiment[order(performancetable$performance)])

rownames(performancetable) = performancetable$experiment

performancetable$experiment = gsub('MNB', 'NB', performancetable$experiment)
performancetable$experiment = gsub('_', ' ', performancetable$experiment)
performancetable$Model = substr(performancetable$experiment, 1, 2)
performancetable$Model = gsub('NB', 'Naive Bayes', performancetable$Model)
performancetable$Model = gsub('NN', 'Neural Network', performancetable$Model)
performancetable$Model = gsub('RF', 'Random Forest', performancetable$Model)

step1_models = rownames(performancetable)[!rownames(performancetable) %in% 
                                            rownames(performancetable)[grep('reduced|qval', rownames(performancetable))]]
step1_models = c(step1_models, 'RF_reducedV3.4_nfeatures21', 'MNB_reducedV3.4_nfeatures21', 'NN_reducedV3.4_nfeatures21')

current_table = performancetable[step1_models, ]

rownames(current_table) = gsub('reducedV...', 'R', rownames(current_table))
rownames(current_table) = gsub('full.features', 'F', rownames(current_table))
rownames(current_table) = gsub('nfeatures', '', rownames(current_table))
rownames(current_table) = gsub('pca\\d+', 'PCA', rownames(current_table))
rownames(current_table) = gsub('_', '.', rownames(current_table))

current_table$experiment = rownames(current_table)
current_table$experiment <- 
  factor(current_table$experiment, levels = current_table$experiment[order(current_table$performance)])

current_colors = c("#00BA38","#E3140F","#497DE7")
output_fn = 'plots/combined_performance_step1'
widths = c(1.5,1,1)

p1 = ggplot(data=current_table, aes(x=experiment, y=accuracyAll, color=Model)) + 
  scale_y_continuous(limit = c((floor(10*min(current_table$performance)) / 10) - 0.1, 1), breaks=seq(0,1,by=.1)) +
  geom_point() +
  theme(plot.title = element_text(hjust = 0.5),
        title = element_text(size=20)) +
  xlab("Model") +
  ylab("Accuracy") +
  theme_bw() + 
  theme(panel.grid = element_line(colour = "#e8e8e8")) +
  geom_errorbar(aes(ymax = current_table$lowerAcc, ymin = current_table$upperAcc), width=0) +
  theme(axis.text.x = element_text(colour="grey20",size=10,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="grey20",size=10,angle=0,vjust=0,face="plain"),  
        axis.title.x = element_text(colour="grey20",size=15,angle=0,hjust=.5,vjust=0,face="plain"),
        axis.title.y = element_blank(),
        legend.position='none') +
  scale_colour_manual(values = current_colors) +
  coord_flip()

p2 = ggplot(data=current_table, aes(x=experiment, y=Kappa, color=Model)) + 
  scale_y_continuous(limit = c((floor(10*min(current_table$performance)) / 10) - 0.1, 1), breaks=seq(0,1,by=.1)) +
  geom_point() +
  theme(plot.title = element_text(hjust = 0.5),
        title = element_text(size=20)) +
  xlab("Model") +
  ylab("K") +
  theme_bw() + 
  theme(panel.grid = element_line(colour = "#e8e8e8")) +
  geom_errorbar(aes(ymax = current_table$upperKappa, ymin = current_table$lowerKappa), width=0) +
  theme(axis.text.x = element_text(colour="grey20",size=10,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_blank(),  
        axis.title.x = element_text(colour="grey20",size=15,angle=0,hjust=.5,vjust=0,face="plain"),
        axis.title.y = element_blank(),
        legend.position = "none") +
  scale_colour_manual(values=current_colors) +
  coord_flip()

p3 = ggplot(data=current_table, aes(x=experiment, y=performance, color=Model)) + 
  scale_y_continuous(limit = c((floor(10*min(current_table$performance)) / 10) - 0.1, 1), breaks=seq(0,1,by=.1)) +
  geom_point() +
  theme(plot.title = element_text(hjust = 0.5),
        title = element_text(size=20)) +
  xlab("Model") +
  ylab("Performance") +
  theme_bw() + 
  theme(panel.grid = element_line(colour = "#e8e8e8")) +
  geom_errorbar(aes(ymax = current_table$upperPerformance, ymin = current_table$lowerPerformance), width=0) +
  theme(axis.text.x = element_text(colour="grey20",size=10,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_blank(),  
        axis.title.x = element_text(colour="grey20",size=15,angle=0,hjust=.5,vjust=0,face="plain"),
        axis.title.y = element_blank(),
        legend.position = c(.26, 0.92),
        legend.title = element_text(size=9)) +
  scale_colour_manual(values=current_colors) +
  geom_hline(yintercept=max(current_table$performance) - 
               (max(current_table$performance) - max(current_table$lowerPerformance)) * 2, linetype="dashed", 
             color = "black", size=0.3) +
  coord_flip()

p4 = plot_grid(p1,p2,p3, nrow=1, rel_widths = widths)
ggsave(paste(output_fn, '.jpeg', sep=''), p4, width=11, height=8)
ggsave(paste(output_fn, '.pdf', sep=''), p4, width=11, height=8)