rm(list = ls())
source("src_R/load_libraries.R")
#library(ggsankey)
library(ggalluvial)

CONFIDENCE_THRESHOLD = 0.8

resultspurity = read.csv('reduce_purity_experiment/resultstable_reducepurity.txt', sep='\t', row.names=1)
labels = read.csv('data_tables/confidence_tables/baseline_probabilities.connectivity_based.sensitivity_power2.Aug_17_2022.tsv', sep='\t', row.names=1)
labels = labels[rownames(resultspurity), ]
countspurity = read.csv('reduce_purity_experiment/countstable_reducepurity.txt', sep='\t')

resultspurity = cbind(labels$confidence, resultspurity)
resultspurity = cbind(labels$cluster, resultspurity)
colnames(resultspurity)[1] = 'true_label'
colnames(resultspurity)[2] = 'target_confidence'
# resultspurity = as_tibble(resultspurity)
# resultspurity_sub = resultspurity[c('initial_pred', 'clus0.6', 'clus0.5', 'clus0.4', 'clus0.3', 'clus0.2', 
#                                     'clus0.15', 'clus0.1', 'clus0.05', 'clus0.02')]

clus_cols =  c('true_label', 'initial_pred', 'clus0.5', 'clus0.4', 'clus0.3', 'clus0.2','clus0.1', 'clus0.02')
conf_cols = c('target_confidence', 'initial_confidence', 'conf0.5', 'conf0.4', 'conf0.3', 'conf0.2','conf0.1', 'conf0.02')

resultspurity_sub = resultspurity[, clus_cols]

confpurity_sub = resultspurity[, conf_cols]

median_confs = apply(confpurity_sub, 2, median)
accuracies = colSums(resultspurity_sub == labels$cluster) / nrow(resultspurity_sub)

resultspurity_sub_lc = resultspurity[resultspurity$initial_confidence <= CONFIDENCE_THRESHOLD, clus_cols]

resultspurity_sub = resultspurity[resultspurity$initial_confidence >= CONFIDENCE_THRESHOLD, clus_cols]

scatter_df = data.frame(median_confs)
colnames(scatter_df) = c('y')
scatter_df$group = 'Median Confidence'
tmp = data.frame(accuracies)
colnames(tmp) = c('y')
tmp$group = 'Accuracy'
scatter_df = rbind(scatter_df, tmp)
scatter_df$purity = rep(colnames(resultspurity_sub), 2)
scatter_df$purity = gsub('clus', 'pur',scatter_df$purity)

#resultspurity_sub_count = resultspurity_sub %>% count(initial_pred, clus0.6, clus0.5, clus0.4, clus0.3, clus0.2, clus0.15, clus0.1, clus0.05, clus0.02)
resultspurity_sub_count = resultspurity_sub %>% count(true_label, initial_pred, clus0.5, 
                                                      clus0.4, clus0.3, clus0.2, clus0.1, clus0.02)

resultspurity_sub_lc_count = resultspurity_sub_lc %>% count(true_label, initial_pred, clus0.5, 
                                                      clus0.4, clus0.3, clus0.2, clus0.1, clus0.02)

colnames(resultspurity_sub_lc_count) = c('Cluster Label', 'Initial Prediction', 'Purity 0.5', 'Purity 0.4', 
                                         'Purity 0.3', 'Purity 0.2', 'Purity 0.1', 'Purity 0.02', 'n')

colnames(resultspurity_sub_count) = c('Cluster Label', 'Initial Prediction', 'Purity 0.5', 'Purity 0.4', 
                                      'Purity 0.3', 'Purity 0.2', 'Purity 0.1', 'Purity 0.02', 'n')

rp_lodes = to_lodes_form(resultspurity_sub_count,
                         key = "Step",
                         axes = 1:(ncol(resultspurity_sub_count) - 1))

rp_lodes_lc = to_lodes_form(resultspurity_sub_lc_count,
                         key = "Step",
                         axes = 1:(ncol(resultspurity_sub_lc_count) - 1))

x_ticks = colnames(resultspurity_sub)

uq_lc = sapply(resultspurity_sub_lc, unique)
uq_lc = sapply(uq_lc, sort, decreasing=T)
stratum_fill_lc = c()
for(col in uq_lc){stratum_fill_lc = c(stratum_fill_lc, col)}
stratum_fill_lc = mapvalues(stratum_fill_lc, from=c(1,2,3,4,5), to=c("#803e98", "#00a4d1", "#f29123", "#4d872d", "#dd2929"))

p = ggplot(data = rp_lodes,
           aes(x = Step, stratum = stratum, alluvium = alluvium,
               y = n, label = stratum)) +
  geom_alluvium(aes(fill = stratum)) +
  geom_stratum(color = "black", size = 0.3,
               fill=rep(rev(c("#803e98", "#00a4d1", "#f29123", "#4d872d", "#dd2929")),
                        (ncol(resultspurity_sub_count) - 1)),
               alpha=1.0) +
  geom_text(stat = "stratum") +
  theme_bw() +
  scale_fill_manual(values = c("#803e98", "#00a4d1", "#f29123", "#4d872d", "#dd2929"), name="Cluster Label") +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5, size=20)
  ) +
  ggtitle('Purity Reduction Simulation - High Confidence Samples')

p_lc = ggplot(data = rp_lodes_lc,
           aes(x = Step, stratum = stratum, alluvium = alluvium,
               y = n, label = stratum)) +
  geom_alluvium(aes(fill = stratum), alpha=0.6) +
  geom_stratum(color = 'black', size = 0.3,
              fill=rep(rev(c("#803e98", "#00a4d1", "#f29123", "#4d872d", "#dd2929")),
                       (ncol(resultspurity_sub_count) - 1)),
               alpha=1.0) +
  geom_text(stat = "stratum") +
  theme_bw() +
  scale_fill_manual(values = c("#803e98", "#00a4d1", "#f29123", "#4d872d", "#dd2929"), name="Cluster Label") +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5, size=20)
  ) +
  ggtitle('Purity Reduction Simulation - Low Confidence Samples')

# p = ggplot(data = rp_lodes,
#            aes(x = Step, stratum = stratum, alluvium = alluvium,
#                y = n, label = stratum)) +
#   geom_alluvium(aes(fill = stratum)) +
#   geom_stratum(color = "white", size = 0.1,
#                alpha=0.4) + 
#   geom_text(stat = "stratum") +
#   theme_bw() +
#   scale_fill_manual(values = c("#803e98", "#00a4d1", "#f29123", "#4d872d", "#dd2929"), name="Initial Cluster") +
#   theme(axis.text.y=element_blank(),
#         axis.ticks.y=element_blank(),
#         axis.title.y = element_blank()
#   )

ggsave('plots/pur_sim_alluvial.png', p, width = 12, height=8)
ggsave('plots/paper_figures/pur_sim_alluvial.pdf', p, width = 12, height=8)

ggsave('plots/pur_lc_sim_alluvial.png', p_lc, width = 12, height=8)
ggsave('plots/paper_figures/pur_lc_sim_alluvial.pdf', p_lc, width = 12, height=8)

level_order = scatter_df$purity[1:8]

melted_df = melt(countspurity)
colnames(melted_df) = c('Event_Type', 'Step', 'Count')

p1 = ggplot(data=melted_df, aes(x=Step, y=Count, group=Event_Type, color=Event_Type)) +
  geom_point() +
  geom_line() +
  theme_bw() +
  theme(panel.grid = element_line(colour = "#C1C1C1")) +
  theme(axis.text.y = element_text(size=10),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=11),
        legend.position = c(0.95, 0.65),
        legend.title = element_text(size = 7), 
        legend.text = element_text(size = 7),
        legend.key.size = unit(3, 'mm'),
        panel.grid.minor.y = element_blank()) +
  scale_color_manual(breaks = countspurity$X,
                     values = c('#000000', '#9c9c9c', 
                                '#039c40', 
                                '#ff0000', '#f77e7e',
                                '#055df5', '#75a8ff')) +
  ggtitle('Accuracy & Confidence vs Purity - High Confidence Samples')

p2 = ggplot(data = scatter_df,
           aes(x=purity, y=y, color=group)) +
  geom_point() +
  theme_bw() +
  scale_x_discrete(limits = level_order) +
  theme(axis.title.y = element_blank(),
        legend.position = c(0.92, 0.86))

p3 = plot_grid(p1, p2, nrow=2, rel_widths = c(1,1), rel_heights = c(1, 3), align='v', axis='l')

ggsave('plots/purity_acc_conf.png', p3, width = 8, height=6)
ggsave('plots/paper_figures/purity_acc_conf.pdf', p3, width=8, height=6)

swapped_samples = rownames(resultspurity)[resultspurity$true_label != resultspurity$initial_pred]
swapped_samples = as_tibble(resultspurity[swapped_samples, ], rownames = NA)
