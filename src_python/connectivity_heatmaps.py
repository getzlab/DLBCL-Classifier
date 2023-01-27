import pandas as pd
from matplotlib import pyplot as plt

VERSION = 'Sep_23_2022'

for k in range(2, 9):
    curr_k = str(k)
    k_conn = pd.read_csv('../clustering_runs/full_dir_combined/k' + curr_k + '/GSM699_cluster_' + VERSION + '.k' + curr_k + '.connectivity.matrix.txt', sep='\t')
    k_labels = pd.read_csv('../clustering_runs/full_dir_combined/k' + curr_k + '/GSM699_cluster_' + VERSION + '.k' + curr_k + '.clustering', sep='\t')
    k_labels = k_labels.sort_values(by=k_labels.columns[-1])

    k_conn = k_conn.loc[k_labels.index, k_labels.index]

    ax = plt.gca()
    im = ax.imshow(k_conn, cmap='gray', interpolation='nearest')
    ax.set_title('k = ' + curr_k)
    cbar = ax.figure.colorbar(im, ax=ax)
    cbar.ax.set_ylabel('connectivity', rotation=-90, va="bottom")
    plt.savefig('../plots/heatmaps/k' + curr_k + '_connectivity_heatmap.pdf')
    plt.savefig('../plots/heatmaps/k' + curr_k + '_connectivity_heatmap.png')
    plt.clf()

for k in range(5, 6):
    curr_k = str(k)
    k_conn = pd.read_csv('../clustering_runs/full_dir_shipp/k' + curr_k + '/GSM699_cluster_' + VERSION + '_SHIPP' + '.k' + curr_k + '.connectivity.matrix.txt',
                         sep='\t')
    k_labels = pd.read_csv('../clustering_runs/full_dir_shipp/k' + curr_k + '/GSM699_cluster_' + VERSION + '_SHIPP' + '.k' + curr_k + '.clustering',
                           sep='\t')
    k_labels = k_labels.sort_values(by=k_labels.columns[-1])

    k_conn = k_conn.loc[k_labels.index, k_labels.index]

    ax = plt.gca()
    im = ax.imshow(k_conn, cmap='gray', interpolation='nearest')
    ax.set_title('k = ' + curr_k)
    cbar = ax.figure.colorbar(im, ax=ax)
    cbar.ax.set_ylabel('connectivity', rotation=-90, va="bottom")
    plt.savefig('../plots/heatmaps/k' + curr_k + '_SHIPP_connectivity_heatmap.pdf')
    plt.savefig('../plots/heatmaps/k' + curr_k + '_SHIPP_connectivity_heatmap.png')
    plt.clf()

for k in range(5, 6):
    curr_k = str(k)
    k_conn = pd.read_csv('../clustering_runs/full_dir_staudt/k' + curr_k + '/GSM699_cluster_' + VERSION + '_STAUDT' + '.k' + curr_k + '.connectivity.matrix.txt',
                         sep='\t')
    k_labels = pd.read_csv('../clustering_runs/full_dir_staudt/k' + curr_k + '/GSM699_cluster_' + VERSION + '_STAUDT' + '.k' + curr_k + '.clustering',
                           sep='\t')
    k_labels = k_labels.sort_values(by=k_labels.columns[-1])

    k_conn = k_conn.loc[k_labels.index, k_labels.index]

    ax = plt.gca()
    im = ax.imshow(k_conn, cmap='gray', interpolation='nearest')
    ax.set_title('k = ' + curr_k)
    cbar = ax.figure.colorbar(im, ax=ax)
    cbar.ax.set_ylabel('connectivity', rotation=-90, va="bottom")
    plt.savefig('../plots/heatmaps/k' + curr_k + '_STAUDT_connectivity_heatmap.pdf')
    plt.savefig('../plots/heatmaps/k' + curr_k + '_STAUDT_connectivity_heatmap.png')
    plt.clf()
