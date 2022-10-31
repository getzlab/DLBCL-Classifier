rm(list = ls())
source("src_R/load_libraries.R")

resultsDropout = read.csv('random_dropout_experiment/resultstable_random_dropout_testset.txt', sep='\t')
resultsAddIn = read.csv('random_add_in_experiment/resultstable_random_add_in_testset.txt', sep='\t')

#####
# Add-in experiment limits: As FP -> 1, TN -> 0, FP/FP = 1. x-axis scale 0->1.
# Alternatively, specificity also mathematically works, but scale must be reversed:
# TP / (TP + FN): As TP -> 0, FN -> 1. 0 / FN = 0. x-axis scale 1->0
#####

#####
# Dropout experiment limits: As TP -> 0, FN -> 1, 0/FN = 0. x-axis scale 1->0.
#####

# FDR = FP / (FP + TN)
# I don't think FDR models our dropout experiment, since FP is just always zero here. We aren't
# introducing any False Positives, only False Negatives.

dropout_label = 'Sensitivity: TP / (TP + FN)'
add_in_label = 'False Positive Rate : FP / (FP + TN)'

newDF = data.frame(resultsDropout$dropout_probability,
                   resultsDropout$mut_count_total / (resultsDropout$mut_count_total + resultsDropout$mut_count_dropped),
                   resultsDropout$amp_count_total / (resultsDropout$amp_count_total + resultsDropout$amp_count_dropped),
                   resultsDropout$del_count_total / (resultsDropout$del_count_total + resultsDropout$del_count_dropped),
                   resultsDropout$sv_count_total / (resultsDropout$sv_count_total + resultsDropout$sv_count_dropped))
colnames(newDF) = c('dropout_probability', 'perc_mut', 'perc_amp', 'perc_del', 'perc_sv')
melted_df = melt(newDF, id = 'dropout_probability')
melted_df$MarkerType = 'Mutation'
melted_df[melted_df$variable == 'perc_amp', 'MarkerType'] = 'Amp'
melted_df[melted_df$variable == 'perc_del', 'MarkerType'] = 'Del'
melted_df[melted_df$variable == 'perc_sv', 'MarkerType'] = 'SV'

resultsDropout_2 = resultsDropout[, c('dropout_probability', 'accuracyAll',	'Kappa', 'Performance',	
                                      'lowerAccuracy',	'upperAccuracy',	'lowerKappa',	'upperKappa',	
                                      'lowerPerformance',	'upperPerformance')]

resultsDropout_2 = rbind(resultsDropout_2)

p = ggplot(resultsDropout_2, aes(x=1-dropout_probability, y=Performance)) +
  geom_point() +
  theme_bw() +
  geom_errorbar(aes(ymax = resultsDropout_2$upperPerformance, ymin = resultsDropout_2$lowerPerformance), width=0.01) +
  labs(x=dropout_label, y='Performance') +
  theme(axis.title.x = element_text(size=22),
        axis.title.y = element_text(size=22),
        axis.text.x = element_text(size=18),
        axis.text.y = element_text(size=18)) +
  ylim(c(0, 1)) +
  ggtitle('Dropout') +
  scale_x_reverse(breaks=seq(0, 1, 0.1), expand = c(0.01, 0), limits=c(1, 0))

p2 = ggplot(data=melted_df, aes(x=dropout_probability, y=value, color=MarkerType)) +
  geom_point() +
  geom_line() +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size=10),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=11),
        legend.position = c(0.2, 0.3),
        legend.title = element_text(size = 7), 
        legend.text = element_text(size = 7),
        legend.key.size = unit(3, 'mm'),
        panel.grid.minor.y = element_blank()) +
  labs(y = 'Relative Saturation') +
  scale_x_continuous(expand = c(0.01, 0), limits = c(0, 1), breaks=seq(0, 1, 0.1)) +
  scale_y_continuous(limits = c(0, 1), breaks=seq(0, 1, 0.2)) +
  scale_color_manual(breaks = c('Mutation', 'Amp', 'Del', 'SV'),
                     values = c('#000000', '#fa2f4e', '#21f4ff', '#00b05e'))

p3 = plot_grid(p2, p, nrow=2, rel_widths = c(1,1), rel_heights = c(1, 3), align='v', axis='l')

ggsave("plots/test_set/performance_dropout_testset.jpeg", p3, width=12, height=8)
ggsave("plots/test_set/performance_dropout_testset.pdf", p3)

p = ggplot(resultsDropout, aes(x=1-dropout_probability, y=meanconfidence)) +
  geom_point() +
  theme_bw() +
  geom_errorbar(aes(ymax = mapply(function(num1) min(c(num1, 1)), resultsDropout$meanconfidence + resultsDropout$stdconfidence), 
                    ymin = resultsDropout$meanconfidence - resultsDropout$stdconfidence), width=0.01) +
  labs(x=dropout_label, y='Mean Confidence') + 
  theme(axis.title.x = element_text(size=22),
        axis.title.y = element_text(size=22),
        axis.text.x = element_text(size=18),
        axis.text.y = element_text(size=18)) +
  ylim(c(0, 1)) +
  scale_x_reverse()

ggsave("plots/test_set/confidence_dropout_testset.jpeg", p)
ggsave("plots/test_set/confidence_dropout_testset.pdf", p)

p = ggplot(resultsDropout, aes(x=1-dropout_probability, y=accuracyAll)) +
  geom_errorbar(aes(ymax = resultsDropout$upperAcc, ymin = resultsDropout$lowerAcc), width=0.01) +
  scale_x_reverse()+
  theme_bw() +
  labs(x=dropout_label, y='Accuracy') + 
  theme(axis.title.x = element_text(size=22),
        axis.title.y = element_text(size=22),
        axis.text.x = element_text(size=18),
        axis.text.y = element_text(size=18)) +
  ylim(c(0, 1)) +
  geom_point()

ggsave("plots/test_set/accuracy_dropout_testset.jpeg", p)
ggsave("plots/test_set/accuracy_dropout_testset.pdf", p)

p = ggplot(resultsDropout, aes(x=1-dropout_probability, y=Kappa)) +
  geom_errorbar(aes(ymax = resultsDropout$upperKappa, ymin = resultsDropout$lowerKappa), width=0.01) +
  labs(x=dropout_label, y='Kappa') +
  scale_x_reverse()+
  theme_bw() +
  ylim(c(0, 1)) +
  theme(axis.title.x = element_text(size=22),
        axis.title.y = element_text(size=22),
        axis.text.x = element_text(size=18),
        axis.text.y = element_text(size=18)) +
  geom_point()


ggsave("plots/test_set/kappa_dropout_testset.jpeg", p)
ggsave("plots/test_set/kappa_dropout_testset.pdf", p)

# Add in experiment 

p = ggplot(resultsAddIn, aes(x=addInFraction, y=Performance)) +
  geom_point() +
  geom_errorbar(aes(ymax = resultsAddIn$upperPerformance, ymin = resultsAddIn$lowerPerformance), width=0.01) +
  labs(x=add_in_label, y='Performance') + 
  ylim(c(0, 1)) +
  theme_bw() +
  theme(axis.title.x = element_text(size=22),
        axis.title.y = element_text(size=22),
        axis.text.x = element_text(size=18),
        axis.text.y = element_text(size=18)) +
  ggtitle('Add-In') +
  scale_x_continuous(expand = c(0.01, 0), limits = c(0, 1), breaks=seq(0, 1, 0.1))

newDF = data.frame(resultsAddIn$addInFraction,
                   resultsAddIn$mut_count_total / (resultsAddIn$mut_count_total - resultsAddIn$mut_count_added),
                   resultsAddIn$amp_count_total / (resultsAddIn$amp_count_total - resultsAddIn$amp_count_added),
                   resultsAddIn$del_count_total / (resultsAddIn$del_count_total - resultsAddIn$del_count_added),
                   resultsAddIn$sv_count_total / (resultsAddIn$sv_count_total - resultsAddIn$sv_count_added))
colnames(newDF) = c('add_in_fraction', 'perc_mut', 'perc_amp', 'perc_del', 'perc_sv')
melted_df = melt(newDF, id = 'add_in_fraction')
melted_df$MarkerType = 'Mutation'
melted_df[melted_df$variable == 'perc_amp', 'MarkerType'] = 'Amp'
melted_df[melted_df$variable == 'perc_del', 'MarkerType'] = 'Del'
melted_df[melted_df$variable == 'perc_sv', 'MarkerType'] = 'SV'

p2 = ggplot(data=melted_df, aes(x=add_in_fraction, y=value, color=MarkerType)) +
  geom_point() +
  geom_line() +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size=10),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=11),
        legend.position = c(0.2, 0.6),
        legend.title = element_text(size = 7), 
        legend.text = element_text(size = 7),
        legend.key.size = unit(3, 'mm'),
        panel.grid.minor.y = element_blank()) +
  labs(y = 'Relative Saturation') +
  scale_x_continuous(expand = c(0.01, 0), limits = c(0, 1), breaks=seq(0, 1, 0.1)) +
  scale_y_continuous(limits = c(1, 12), breaks=seq(1,12,2)) +
  scale_color_manual(breaks = c('Mutation', 'Amp', 'Del', 'SV'),
                     values = c('#000000', '#fa2f4e', '#21f4ff', '#00b05e'))

p3 = plot_grid(p2, p, nrow=2, rel_widths = c(1,1), rel_heights = c(1, 3), align='v', axis='l')

ggsave("plots/test_set/performance_addIn_testset.jpeg", p3, width = 12, height=8)
ggsave("plots/test_set/performance_addIn_testset.pdf", p3)

p = ggplot(resultsAddIn, aes(x=addInFraction, y=accuracyAll)) +
  geom_errorbar(aes(ymax = resultsAddIn$upperAcc, ymin = resultsAddIn$lowerAcc), width=0.01) +
  geom_point() +
  labs(x=add_in_label, y='Accuracy') + 
  ylim(c(0, 1)) +
  theme_bw() +
  theme(axis.title.x = element_text(size=22),
        axis.title.y = element_text(size=22),
        axis.text.x = element_text(size=18),
        axis.text.y = element_text(size=18))


ggsave("plots/test_set/accuracy_addIn_testset.jpeg", p)
ggsave("plots/test_set/accuracy_addIn_testset.pdf", p)

p = ggplot(resultsAddIn, aes(x=addInFraction, y=Kappa)) +
  ylim(c(0, 1)) +
  theme_bw() +
  labs(x=add_in_label, y='Kappa') + 
  theme(axis.title.x = element_text(size=22),
        axis.title.y = element_text(size=22),
        axis.text.x = element_text(size=18),
        axis.text.y = element_text(size=18)) +
  geom_errorbar(aes(ymax = resultsAddIn$upperKappa, ymin = resultsAddIn$lowerKappa), width=0.01) +
  geom_point()

ggsave("plots/test_set/kappa_addIn_testset.jpeg", p)
ggsave("plots/test_set/kappa_addIn_testset.pdf", p)


p = ggplot(resultsAddIn, aes(x=addInFraction, y=meanconfidence)) +
  geom_point() +
  theme_bw() +
  geom_errorbar(aes(ymax = mapply(function(num1) min(c(num1, 1)), resultsAddIn$meanconfidence + resultsAddIn$stdconfidence), 
                    ymin = resultsAddIn$meanconfidence - resultsAddIn$stdconfidence), width=0.01) +
  labs(x=add_in_label, y='Mean Confidence') + 
  theme(axis.title.x = element_text(size=22),
        axis.title.y = element_text(size=22),
        axis.text.x = element_text(size=18),
        axis.text.y = element_text(size=18)) +
  ylim(c(0, 1))

ggsave("plots/test_set/confidence_addIn_testset.jpeg", p)
ggsave("plots/test_set/confidence_addIn_testset.pdf", p)


