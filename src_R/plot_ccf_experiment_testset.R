rm(list = ls())
source("src_R/load_libraries.R")

resultsCCF = read.csv('ccf_threshold_experiment/resultstable_ccf_threshold_test.txt', sep='\t')

ccf_label = 'CCF Threshold (>=)'

p_perf = ggplot(resultsCCF, aes(x=ccf_threshold, y=Performance)) +
  geom_point() +
  theme_bw() +
  geom_errorbar(aes(ymax = resultsCCF$upperPerformance, ymin = resultsCCF$lowerPerformance), width=0.01) +
  labs(x=ccf_label, y='Performance') + 
  theme(axis.title.x = element_text(size=22),
        axis.title.y = element_text(size=22),
        axis.text.x = element_text(size=18),
        axis.text.y = element_text(size=18)) +
  scale_x_continuous(expand = c(0.01, 0), breaks = seq(0, 1, 0.10), limits=c(0,1)) +
  scale_y_continuous(expand = c(0.01, 0), breaks = seq(0, 1, 0.10), limits=c(0,1))

p_conf = ggplot(resultsCCF, aes(x=ccf_threshold, y=meanconfidence)) +
  geom_point() +
  theme_bw() +
  geom_errorbar(aes(ymax = resultsCCF$meanconfidence + resultsCCF$stdconfidence, 
                    ymin = resultsCCF$meanconfidence - resultsCCF$stdconfidence), width=0.01) +
  labs(x=ccf_label, y='Mean Confidence') + 
  theme(axis.title.x = element_text(size=22),
        axis.title.y = element_text(size=22),
        axis.text.x = element_text(size=18),
        axis.text.y = element_text(size=18)) +
  ylim(c(0, 1))

p_acc = ggplot(resultsCCF, aes(x=ccf_threshold, y=accuracyAll)) +
  geom_errorbar(aes(ymax = resultsCCF$upperAcc, ymin = resultsCCF$lowerAcc), width=0.01) +
  theme_bw() +
  labs(x=ccf_label, y='Accuracy') + 
  theme(axis.title.x = element_text(size=22),
        axis.title.y = element_text(size=22),
        axis.text.x = element_text(size=18),
        axis.text.y = element_text(size=18)) +
  ylim(c(0, 1)) +
  geom_point()

p_calib = ggplot(resultsCCF, aes(x=ccf_threshold, y=Kappa)) +
  geom_errorbar(aes(ymax = resultsCCF$upperKappa, ymin = resultsCCF$lowerKappa), width=0.01) +
  labs(x=ccf_label, y='Kappa') +
  theme_bw() +
  ylim(c(0, 1)) +
  theme(axis.title.x = element_text(size=22),
        axis.title.y = element_text(size=22),
        axis.text.x = element_text(size=18),
        axis.text.y = element_text(size=18)) +
  geom_point()


newDF = data.frame(resultsCCF$ccf_threshold,
                   resultsCCF$mut_count_kept / (resultsCCF$mut_count_dropped + resultsCCF$mut_count_kept),
                   resultsCCF$amp_count_kept / (resultsCCF$amp_count_dropped + resultsCCF$amp_count_kept),
                   resultsCCF$del_count_kept / (resultsCCF$del_count_dropped + resultsCCF$del_count_kept),
                   resultsCCF$sv_count_kept / (resultsCCF$sv_count_dropped + resultsCCF$sv_count_kept))
colnames(newDF) = c('ccf_threshold', 'perc_mut', 'perc_amp', 'perc_del', 'perc_sv')
melted_df = melt(newDF, id = 'ccf_threshold')
melted_df$MarkerType = 'Mutation'
melted_df[melted_df$variable == 'perc_amp', 'MarkerType'] = 'Amp'
melted_df[melted_df$variable == 'perc_del', 'MarkerType'] = 'Del'
melted_df[melted_df$variable == 'perc_sv', 'MarkerType'] = 'SV'

p2 = ggplot(data=melted_df, aes(x=ccf_threshold, y=value, color=MarkerType)) +
  geom_point() +
  geom_line() +
  theme_bw() +
  theme(panel.grid = element_line(colour = "#C1C1C1")) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size=10),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=11),
        legend.position = c(0.3, 0.3),
        legend.title = element_text(size = 7), 
        legend.text = element_text(size = 7),
        legend.key.size = unit(3, 'mm'),
        panel.grid.minor.y = element_blank()) +
  labs(y = 'Fraction Remaining') +
  scale_x_continuous(expand = c(0.01, 0), limits = c(0, 1), breaks=seq(0, 1, 0.1)) +
  scale_y_continuous(expand = c(0.01, 0), breaks = seq(0, 1, 0.10), limits=c(0,1)) +
  scale_color_manual(breaks = c('Mutation', 'Amp', 'Del', 'SV'),
                     values = c('#000000', '#fa2f4e', '#21f4ff', '#00b05e'))

p3 = plot_grid(p2, p_perf, nrow=2, rel_widths = c(1,1), rel_heights = c(1, 3), align='v', axis='l')

ggsave("plots/test_set/kappa_ccf_testset.jpeg", p_calib)
ggsave("plots/test_set/kappa_ccf_testset.pdf", p_calib)
ggsave("plots/test_set/accuracy_ccf_testset.jpeg", p_acc)
ggsave("plots/test_set/accuracy_ccf_testset.pdf", p_acc)
ggsave("plots/test_set/performance_ccf_testset.jpeg", p3)
ggsave("plots/test_set/performance_ccf_testset.pdf", p3, width=10, height=8)
ggsave("plots/test_set/confidence_ccf_testset.jpeg", p_conf)
ggsave("plots/test_set/confidence_ccf_testset.pdf", p_conf)

