rm(list = ls())
source('src_R/load_libraries.R')

performancetable = read.csv('evaluation_panel_sets/evaluation_panels.tsv', sep='\t', row.names=1)

current_table = performancetable

current_table$experiment = rownames(current_table)
current_table$experiment <- 
  factor(current_table$experiment, levels = current_table$experiment[order(current_table$performance)])

current_table$Model = 'Panel'
current_table['dlbclass', 'Model'] = 'DLBclass'

current_colors = c("#fc0303","#000000")
output_fn = 'plots/combined_performance_panels'
widths = c(1.5,1,1)

p1 = ggplot(data=current_table, aes(x=experiment, y=accuracy, color=Model)) + 
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

p2 = ggplot(data=current_table, aes(x=experiment, y=kappa, color=Model)) + 
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
  coord_flip()

p4 = plot_grid(p1,p2,p3, nrow=1, rel_widths = widths)
ggsave(paste(output_fn, '.jpeg', sep=''), p4, width=8, height=2)
ggsave(paste(output_fn, '.pdf', sep=''), p4, width=8, height=2)