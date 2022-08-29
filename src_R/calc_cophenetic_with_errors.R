rm(list = ls())
source('src_R/load_libraries.R')

NUM_ITER = 200

conn_mat_2 = read.csv('clustering_runs/full_dir_combined/k2/GSM699_cluster_Aug_17_2022.k2.connectivity.matrix.txt', sep='\t')
conn_mat_3 = read.csv('clustering_runs/full_dir_combined/k3/GSM699_cluster_Aug_17_2022.k3.connectivity.matrix.txt', sep='\t')
conn_mat_4 = read.csv('clustering_runs/full_dir_combined/k4/GSM699_cluster_Aug_17_2022.k4.connectivity.matrix.txt', sep='\t')
conn_mat_5 = read.csv('clustering_runs/full_dir_combined/k5/GSM699_cluster_Aug_17_2022.k5.connectivity.matrix.txt', sep='\t')
conn_mat_6 = read.csv('clustering_runs/full_dir_combined/k6/GSM699_cluster_Aug_17_2022.k6.connectivity.matrix.txt', sep='\t')
conn_mat_7 = read.csv('clustering_runs/full_dir_combined/k7/GSM699_cluster_Aug_17_2022.k7.connectivity.matrix.txt', sep='\t')
conn_mat_8 = read.csv('clustering_runs/full_dir_combined/k8/GSM699_cluster_Aug_17_2022.k8.connectivity.matrix.txt', sep='\t')

c_mats = list(conn_mat_2, conn_mat_3, conn_mat_4, conn_mat_5, conn_mat_6, conn_mat_7, conn_mat_8)

all_rho_vals = c()
rho_means = c()
rho_sds = c()

for (cm in c_mats){
  all_rho = c()
  for (i in 1:NUM_ITER){
    resamp = sample(1:nrow(cm), nrow(cm), replace=TRUE)
    conn_mat_bstrap <- cm[resamp, resamp]
    
    dis_mat = as.dist(1-conn_mat_bstrap)
    
    HC <- hclust(dis_mat, method="average")
    dis_coph <- cophenetic(HC)
    rho <- signif((cor(dis_mat, dis_coph)), digits = 4)
    all_rho = c(all_rho, rho)
  }
  
  rho_means = c(rho_means, mean(all_rho))
  rho_sds = c(rho_sds, sd(all_rho))
  all_rho_vals = c(all_rho_vals, all_rho)
}

current_table = data.frame(seq(2, 8, 1), rho_means, rho_sds)
colnames(current_table)[1] = 'k'
current_table$color = c('black', 'black', 'black', 'red', 'black', 'black', 'black')

p = ggplot(data=current_table, aes(x=k, y=rho_means, color=color)) + 
  geom_point() +
  theme(plot.title = element_text(hjust = 0.5),
        title = element_text(size=20)) +
  xlab("Number of Clusters") +
  ylab("Rho") +
  theme_bw() + 
  ggtitle('Combined') +
  scale_x_discrete(limits=c(2,3,4,5,6,7,8)) +
  ylim(c(0.7, 1)) +
  scale_color_manual(values=c('black', 'red')) +
  theme(panel.grid = element_line(colour = "#e8e8e8")) +
  geom_errorbar(aes(ymax = rho_means + rho_sds, 
                    ymin = rho_means - rho_sds), width=0) +
  theme(axis.text.x = element_text(colour="grey20",size=10,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="grey20",size=10,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="grey20",size=15,angle=0,hjust=.5,vjust=0,face="plain"),
        axis.title.y = element_text(colour="grey20",size=15,angle=0,hjust=.5,vjust=0.5,face="plain"),
        legend.position = "none")

conn_mat_2 = read.csv('clustering_runs/full_dir_staudt/k2/GSM699_cluster_Aug_17_2022_STAUDT.k2.connectivity.matrix.txt', sep='\t')
conn_mat_3 = read.csv('clustering_runs/full_dir_staudt/k3/GSM699_cluster_Aug_17_2022_STAUDT.k3.connectivity.matrix.txt', sep='\t')
conn_mat_4 = read.csv('clustering_runs/full_dir_staudt/k4/GSM699_cluster_Aug_17_2022_STAUDT.k4.connectivity.matrix.txt', sep='\t')
conn_mat_5 = read.csv('clustering_runs/full_dir_staudt/k5/GSM699_cluster_Aug_17_2022_STAUDT.k5.connectivity.matrix.txt', sep='\t')
conn_mat_6 = read.csv('clustering_runs/full_dir_staudt/k6/GSM699_cluster_Aug_17_2022_STAUDT.k6.connectivity.matrix.txt', sep='\t')
conn_mat_7 = read.csv('clustering_runs/full_dir_staudt/k7/GSM699_cluster_Aug_17_2022_STAUDT.k7.connectivity.matrix.txt', sep='\t')
conn_mat_8 = read.csv('clustering_runs/full_dir_staudt/k8/GSM699_cluster_Aug_17_2022_STAUDT.k8.connectivity.matrix.txt', sep='\t')

c_mats = list(conn_mat_2, conn_mat_3, conn_mat_4, conn_mat_5, conn_mat_6, conn_mat_7, conn_mat_8)

all_rho_vals = c()
rho_means = c()
rho_sds = c()

for (cm in c_mats){
  all_rho = c()
  for (i in 1:NUM_ITER){
    resamp = sample(1:nrow(cm), nrow(cm), replace=TRUE)
    conn_mat_bstrap <- cm[resamp, resamp]
    
    dis_mat = as.dist(1-conn_mat_bstrap)
    
    HC <- hclust(dis_mat, method="average")
    dis_coph <- cophenetic(HC)
    rho <- signif((cor(dis_mat, dis_coph)), digits = 4)
    all_rho = c(all_rho, rho)
  }
  
  rho_means = c(rho_means, mean(all_rho))
  rho_sds = c(rho_sds, sd(all_rho))
  all_rho_vals = c(all_rho_vals, all_rho)
}

current_table_staudt = data.frame(seq(2, 8, 1), rho_means, rho_sds)
colnames(current_table_staudt)[1] = 'k'

p_st = ggplot(data=current_table_staudt, aes(x=k, y=rho_means)) + 
  geom_point() +
  theme(plot.title = element_text(hjust = 0.5),
        title = element_text(size=20)) +
  xlab("Number of Clusters") +
  ylab("Rho") +
  theme_bw() + 
  scale_x_discrete(limits=c(2,3,4,5,6,7,8)) +
  ggtitle('Staudt') +
  theme(panel.grid = element_line(colour = "#e8e8e8")) +
  ylim(c(0.7, 1)) +
  geom_errorbar(aes(ymax = current_table_staudt$rho_means + current_table_staudt$rho_sds, 
                    ymin = current_table_staudt$rho_means - current_table_staudt$rho_sds), width=0) +
  theme(axis.text.x = element_text(colour="grey20",size=10,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="grey20",size=10,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="grey20",size=15,angle=0,hjust=.5,vjust=0,face="plain"),
        axis.title.y = element_text(colour="grey20",size=15,angle=0,hjust=.5,vjust=0.5,face="plain"),
        legend.position = "none")


conn_mat_2 = read.csv('clustering_runs/full_dir_shipp/k2/GSM699_cluster_Aug_17_2022_SHIPP.k2.connectivity.matrix.txt', sep='\t')
conn_mat_3 = read.csv('clustering_runs/full_dir_shipp/k3/GSM699_cluster_Aug_17_2022_SHIPP.k3.connectivity.matrix.txt', sep='\t')
conn_mat_4 = read.csv('clustering_runs/full_dir_shipp/k4/GSM699_cluster_Aug_17_2022_SHIPP.k4.connectivity.matrix.txt', sep='\t')
conn_mat_5 = read.csv('clustering_runs/full_dir_shipp/k5/GSM699_cluster_Aug_17_2022_SHIPP.k5.connectivity.matrix.txt', sep='\t')
conn_mat_6 = read.csv('clustering_runs/full_dir_shipp/k6/GSM699_cluster_Aug_17_2022_SHIPP.k6.connectivity.matrix.txt', sep='\t')
conn_mat_7 = read.csv('clustering_runs/full_dir_shipp/k7/GSM699_cluster_Aug_17_2022_SHIPP.k7.connectivity.matrix.txt', sep='\t')
conn_mat_8 = read.csv('clustering_runs/full_dir_shipp/k8/GSM699_cluster_Aug_17_2022_SHIPP.k8.connectivity.matrix.txt', sep='\t')

c_mats = list(conn_mat_2, conn_mat_3, conn_mat_4, conn_mat_5, conn_mat_6, conn_mat_7, conn_mat_8)

all_rho_vals = c()
rho_means = c()
rho_sds = c()

for (cm in c_mats){
  all_rho = c()
  for (i in 1:NUM_ITER){
    resamp = sample(1:nrow(cm), nrow(cm), replace=TRUE)
    conn_mat_bstrap <- cm[resamp, resamp]
    
    dis_mat = as.dist(1-conn_mat_bstrap)
    
    HC <- hclust(dis_mat, method="average")
    dis_coph <- cophenetic(HC)
    rho <- signif((cor(dis_mat, dis_coph)), digits = 4)
    all_rho = c(all_rho, rho)
  }
  
  rho_means = c(rho_means, mean(all_rho))
  rho_sds = c(rho_sds, sd(all_rho))
  all_rho_vals = c(all_rho_vals, all_rho)
}

current_table_shipp = data.frame(seq(2, 8, 1), rho_means, rho_sds)
colnames(current_table_shipp)[1] = 'k'

p_sh = ggplot(data=current_table_shipp, aes(x=k, y=rho_means)) + 
  geom_point() +
  theme(plot.title = element_text(hjust = 0.5),
        title = element_text(size=20)) +
  xlab("Number of Clusters") +
  ylab("Rho") +
  ggtitle("Shipp") +
  theme_bw() + 
  scale_x_discrete(limits=c(2,3,4,5,6,7,8)) +
  theme(panel.grid = element_line(colour = "#e8e8e8")) +
  ylim(c(0.7, 1)) +
  geom_errorbar(aes(ymax = current_table_shipp$rho_means + current_table_shipp$rho_sds, 
                    ymin = current_table_shipp$rho_means - current_table_shipp$rho_sds), width=0) +
  theme(axis.text.x = element_text(colour="grey20",size=10,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="grey20",size=10,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="grey20",size=15,angle=0,hjust=.5,vjust=0,face="plain"),
        axis.title.y = element_text(colour="grey20",size=15,angle=0,hjust=.5,vjust=0.5,face="plain"),
        legend.position = "none")

ggsave('plots/cophenetic_combined.png', p, height=6, width=6)
ggsave('plots/paper_figures/cophenetic_combined.pdf', p, height=6, width=6)
ggsave('plots/cophenetic_shipp.png', p_sh, height=6, width=6)
ggsave('plots/paper_figures/cophenetic_shipp.pdf', p_sh, height=6, width=6)
ggsave('plots/cophenetic_staudt.png', p_st, height=6, width=6)
ggsave('plots/paper_figures/cophenetic_staudt.pdf', p_st, height=6, width=6)

write.table(current_table, 'clustering_runs/cophenetic_with_errors_combined.tsv', sep='\t', row.names = FALSE)
write.table(current_table_staudt, 'clustering_runs/cophenetic_with_errors_STAUDT.tsv', sep='\t', row.names = FALSE)
write.table(current_table_shipp, 'clustering_runs/cophenetic_with_errors_SHIPP.tsv', sep='\t', row.names = FALSE)
