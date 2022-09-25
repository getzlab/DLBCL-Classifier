rm(list = ls())
source("src_R/load_libraries.R")

nm_gsm_features = read.csv('data_tables/NM_supp_cnmf/c)ConsClust. - Marker Selection-Table 1.tsv', sep='\t', header = 2, skip=1)

nm_c1_sig = nm_gsm_features[(nm_gsm_features$Cluster == 1) & (nm_gsm_features$fisher_q_withinCluster < 0.10), 'Alteration']
nm_c2_sig = nm_gsm_features[(nm_gsm_features$Cluster == 2) & (nm_gsm_features$fisher_q_withinCluster < 0.10), 'Alteration']
nm_c3_sig = nm_gsm_features[(nm_gsm_features$Cluster == 3) & (nm_gsm_features$fisher_q_withinCluster < 0.10), 'Alteration']
nm_c4_sig = nm_gsm_features[(nm_gsm_features$Cluster == 4) & (nm_gsm_features$fisher_q_withinCluster < 0.10), 'Alteration']
nm_c5_sig = nm_gsm_features[(nm_gsm_features$Cluster == 5) & (nm_gsm_features$fisher_q_withinCluster < 0.10), 'Alteration']

nm_c1_sig = toupper(make.names(nm_c1_sig))
nm_c2_sig = toupper(make.names(nm_c2_sig))
nm_c3_sig = toupper(make.names(nm_c3_sig))
nm_c4_sig = toupper(make.names(nm_c4_sig))
nm_c5_sig = toupper(make.names(nm_c5_sig))

new_qval_df =  read.csv('data_tables/qval_dfs/fisher_exact_5x2_17-Aug-2022.combined.tsv',
               sep='\t')

new_c1_sig = new_qval_df[((new_qval_df$C1_nf >= 0.30) & (new_qval_df$q < 0.10)) | 
                           ((new_qval_df$q < 0.10) & (new_qval_df$cluster == 'C1')), 'gene']
new_c2_sig = new_qval_df[((new_qval_df$C2_nf >= 0.30) & (new_qval_df$q < 0.10)) | 
                           ((new_qval_df$q < 0.10) & (new_qval_df$cluster == 'C2')), 'gene']
new_c3_sig = new_qval_df[((new_qval_df$C3_nf >= 0.30) & (new_qval_df$q < 0.10)) | 
                           ((new_qval_df$q < 0.10) & (new_qval_df$cluster == 'C3')), 'gene']
new_c4_sig = new_qval_df[((new_qval_df$C4_nf >= 0.30) & (new_qval_df$q < 0.10)) | 
                           ((new_qval_df$q < 0.10) & (new_qval_df$cluster == 'C4')), 'gene']
new_c5_sig = new_qval_df[((new_qval_df$C5_nf >= 0.30) & (new_qval_df$q < 0.10)) | 
                           ((new_qval_df$q < 0.10) & (new_qval_df$cluster == 'C5')), 'gene']

# Manually add in a couple edge cases/exceptions
new_c2_sig = c(new_c2_sig, 'X10Q23.31.DEL')
new_c3_sig = c(new_c3_sig, 'X12P.AMP')
new_c5_sig = c(new_c5_sig, 'ETS1')

pdf("./plots/paper_figures/c1_venn.pdf", height = 11, width = 8.5, paper = "letter")
grid.newpage()
v = venn.diagram(
  x = list(nm_c1_sig, new_c1_sig),
  category.names = c("C1 NM" , "C1 New"),
  scaled=FALSE,
  disable.logging = TRUE,
  filename=NULL,
  cat.pos = c(60, -60),
  cat.dist = c(0.35, 0.35),
  margin=0.1
)

c1_new_only = setdiff(new_c1_sig, nm_c1_sig)
c1_nm_only = setdiff(nm_c1_sig, new_c1_sig)
c1_both = intersect(nm_c1_sig, new_c1_sig)

grid.draw(v)
grid.text(paste(c1_nm_only, collapse='\n'), x = 0.89,  y=0.6,
          gp=gpar(fontsize=10))
grid.text(paste(c1_new_only, collapse='\n'), x = 0.11,  y=0.42,
          gp=gpar(fontsize=10))
grid.text(paste('both\n', paste(c1_both, collapse='\n'), collapse='\n'), x = 0.5,  y=0.18,
          gp=gpar(fontsize=10))

dev.off()

pdf("./plots/paper_figures/c2_venn.pdf", height = 11, width = 8.5, paper = "letter")
grid.newpage()
v = venn.diagram(
  x = list(nm_c2_sig, new_c2_sig),
  category.names = c("C2 NM" , "C2 New"),
  scaled=FALSE,
  disable.logging = TRUE,
  filename=NULL,
  cat.pos = c(60, -60),
  cat.dist = c(0.35, 0.35),
  margin=0.1
)

c2_new_only = setdiff(new_c2_sig, nm_c2_sig)
c2_nm_only = setdiff(nm_c2_sig, new_c2_sig)
c2_both = intersect(nm_c2_sig, new_c2_sig)

grid.draw(v)
grid.text(paste(c2_nm_only, collapse='\n'), x = 0.89,  y=0.6,
          gp=gpar(fontsize=10))
grid.text(paste(c2_new_only, collapse='\n'), x = 0.11,  y=0.42,
          gp=gpar(fontsize=10))
grid.text('both', x=0.5, y=0.3)
grid.text(paste(c2_both[0:as.integer(length(c2_both) / 2)], collapse='\n'), x = 0.4,  y=0.15,
          gp=gpar(fontsize=8))
grid.text(paste(c2_both[as.integer(length(c2_both) / 2):length(c2_both)], collapse='\n'), x = 0.6,  y=0.15,
          gp=gpar(fontsize=8))

dev.off()

grid.newpage()
pdf("./plots/paper_figures/c3_venn.pdf", height = 11, width = 8.5, paper = "letter")
v = venn.diagram(
  x = list(nm_c3_sig, new_c3_sig),
  category.names = c("C3 NM" , "C3 New"),
  scaled=FALSE,
  disable.logging = TRUE,
  filename=NULL,
  cat.pos = c(60, -60),
  cat.dist = c(0.35, 0.35),
  margin=0.1
)

c3_new_only = setdiff(new_c3_sig, nm_c3_sig)
c3_nm_only = setdiff(nm_c3_sig, new_c3_sig)
c3_both = intersect(nm_c3_sig, new_c3_sig)

grid.draw(v)
grid.text(paste(c3_nm_only, collapse='\n'), x = 0.89,  y=0.6,
          gp=gpar(fontsize=10))
grid.text(paste(c3_new_only, collapse='\n'), x = 0.11,  y=0.42,
          gp=gpar(fontsize=10))
grid.text('both', x=0.5, y=0.3)
grid.text(paste(c3_both[0:as.integer(length(c3_both) / 2)], collapse='\n'), x = 0.4,  y=0.15,
          gp=gpar(fontsize=8))
grid.text(paste(c3_both[as.integer(length(c3_both) / 2):length(c3_both)], collapse='\n'), x = 0.6,  y=0.15,
          gp=gpar(fontsize=8))

dev.off()

grid.newpage()
pdf("./plots/paper_figures/c4_venn.pdf", height = 11, width = 8.5, paper = "letter")
v = venn.diagram(
  x = list(nm_c4_sig, new_c4_sig),
  category.names = c("C4 NM" , "C4 New"),
  scaled=FALSE,
  disable.logging = TRUE,
  filename=NULL,
  cat.pos = c(60, -60),
  cat.dist = c(0.35, 0.35),
  margin=0.1
)

c4_new_only = setdiff(new_c4_sig, nm_c4_sig)
c4_nm_only = setdiff(nm_c4_sig, new_c4_sig)
c4_nm_only = gsub('HIST1H2BK', 'HIST1H2BK (Artifact)', c4_nm_only)
c4_both = intersect(nm_c4_sig, new_c4_sig)

grid.draw(v)
grid.text(paste(c4_nm_only, collapse='\n'), x = 0.89,  y=0.6,
          gp=gpar(fontsize=10))
grid.text(paste(c4_new_only, collapse='\n'), x = 0.11,  y=0.42,
          gp=gpar(fontsize=10))
grid.text('both', x=0.5, y=0.3)
grid.text(paste(c4_both[0:as.integer(length(c4_both) / 2)], collapse='\n'), x = 0.4,  y=0.15,
          gp=gpar(fontsize=8))
grid.text(paste(c4_both[as.integer(length(c4_both) / 2):length(c4_both)], collapse='\n'), x = 0.6,  y=0.15,
          gp=gpar(fontsize=8))

dev.off()


grid.newpage()
pdf("./plots/paper_figures/c5_venn.pdf", height = 11, width = 8.5, paper = "letter")
v = venn.diagram(
  x = list(nm_c5_sig, new_c5_sig),
  category.names = c("C5 NM" , "C5 New"),
  scaled=FALSE,
  disable.logging = TRUE,
  filename=NULL,
  cat.pos = c(60, -60),
  cat.dist = c(0.35, 0.35),
  margin=0.1
)

c5_new_only = setdiff(new_c5_sig, nm_c5_sig)
c5_nm_only = setdiff(nm_c5_sig, new_c5_sig)
c5_both = intersect(nm_c5_sig, new_c5_sig)

grid.draw(v)
grid.text(paste(c5_nm_only, collapse='\n'), x = 0.89,  y=0.6,
          gp=gpar(fontsize=10))
grid.text(paste(c5_new_only, collapse='\n'), x = 0.11,  y=0.38,
          gp=gpar(fontsize=10))
grid.text('both', x=0.5, y=0.3)
grid.text(paste(c5_both[0:as.integer(length(c5_both) / 2)], collapse='\n'), x = 0.4,  y=0.15,
          gp=gpar(fontsize=8))
grid.text(paste(c5_both[as.integer(length(c5_both) / 2):length(c5_both)], collapse='\n'), x = 0.6,  y=0.15,
          gp=gpar(fontsize=8))

dev.off()

