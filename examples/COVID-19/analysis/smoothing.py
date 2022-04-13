import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pylab as plt

import skfda
import skfda.misc.hat_matrix as hm
import skfda.preprocessing.smoothing.validation as val
from skfda.preprocessing.smoothing import KernelSmoother

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'Arial'

# create default cmap
def create_cmap(cmap_name, cmap_list):
    cmap_o = plt.get_cmap(cmap_name)
    for i in range(cmap_o.N):
        rgba = cmap_o(i)
        cmap_list.append(mpl.colors.rgb2hex(rgba))
    return cmap_list

cmap = []
for name in ['Paired', 'tab20b', 'tab20c']:
    cmap = create_cmap(name, cmap)

def build_fda(df_zscores, df_metadata, cell_type = 'cMono', perturbation = 'COVID_SEV'):
    # preparation
    df_metadata['time_v2'] = [int(i.split('d')[1]) for i in df_metadata['time']]
    df_metadata = df_metadata.loc[(df_metadata['cell_type'] == cell_type) & (df_metadata['perturbation'] == perturbation), :].copy()
    df_metadata = df_metadata.sort_values(['time_v2'], ascending = True)
    timepoints = list(df_metadata['time_v2'])

    comps = list(df_metadata.index.values)
    df_zscores = df_zscores[comps]

    data_matrix = []
    for gene in df_zscores.index.values:
        gene_zvalues = df_zscores.loc[gene, :]
        data_matrix.append(gene_zvalues)

    # create representation of functional data
    fd = skfda.FDataGrid(
        data_matrix = data_matrix,
        grid_points = timepoints
    )
    genes = df_zscores.index.values
    return fd, genes


def perform_smoothing(fd_whole, df_cluster, genes, n_neighbors = 2, bandwidth = 1, 
                        cluster_key = 'clusters_fuzzy', output_folder = 'figures/smoothing_percluster/'):

    clusters = df_cluster[cluster_key].unique()
    fd_mean = fd_whole.mean() # why it is empty?
    
    print(fd_whole.size)
    print(fd_whole)
    print(fd_mean.size)
    print(fd_mean)
    
    fig, ax = plt.subplots()
    fd_mean.plot(axes = ax, group = [0], group_colors = ['black'], linewidth = 5)
    fig.savefig(output_folder + 'avg_curve.png', bbox_inches = 'tight')   

    for color, cluster in zip(cmap, clusters):
        # initial
        df_subcluster = df_cluster.loc[df_cluster[cluster_key] == cluster, :].copy()
        genes_sub = df_subcluster['genes']
        sub_indices = list(np.where(np.in1d(genes, genes_sub))[0])
        fd = fd_whole[sub_indices]

        # Basic visualization
        fig, axes = plt.subplots(1,2,figsize = (8,3))
        fd.plot(axes = axes[0])
        fd.scatter(s = 0.5, axes = axes[1])
        fig.savefig(output_folder + 'raw_' + str(cluster) + '.pdf', bbox_inches = 'tight')

        # do smoothing
        # K-nearest neighbours kernel smoothing.
        knn = KernelSmoother(kernel_estimator=hm.KNeighborsHatMatrix(n_neighbors = n_neighbors))
        knn.fit(fd)
        knn_fd = knn.transform(fd)

        # Local linear regression kernel smoothing.
        llr = KernelSmoother(kernel_estimator=hm.LocalLinearRegressionHatMatrix(bandwidth = bandwidth))
        llr.fit(fd)
        llr_fd = llr.transform(fd)

        # # Nadaraya-Watson kernel smoothing.
        nw = KernelSmoother(kernel_estimator=hm.NadarayaWatsonHatMatrix(bandwidth = bandwidth))
        nw.fit(fd)
        nw_fd = nw.transform(fd)

        # visualize fitted curves
        fig, ax = plt.subplots()
        knn_fd.plot(axes = ax, group = [0] * fd.size, group_colors = [color], alpha = 0.5)
        fd_mean.plot(axes = ax, group = [0], group_colors = ['black'])
        fig.savefig(output_folder + 'knn_smoothing_' + str(cluster) + '.pdf', bbox_inches = 'tight')    

        fig, ax = plt.subplots()
        llr_fd.plot(axes = ax, group = [0] * fd.size, group_colors = [color], alpha = 0.5)
        fd_mean.plot(axes = ax, group = [0], group_colors = ['black'])
        fig.savefig(output_folder + 'LR_smoothing_' + str(cluster) + '.pdf', bbox_inches = 'tight')    

        fig, ax = plt.subplots()
        nw_fd.plot(axes = ax, group = [0] * fd.size, group_colors = [color], alpha = 0.5)
        fd_mean.plot(axes = ax, group = [0], group_colors = ['black'])
        fig.savefig(output_folder + 'NW_smoothing_' + str(cluster) + '.pdf', bbox_inches = 'tight')    

if __name__ == '__main__':
    # initialization
    df_zscores = pd.read_csv('pairwise_zscores_HVG1000_combined_rep3.txt', sep = '\t', header = 0, index_col = 0)
    df_zscores.values[np.isnan(df_zscores.values)] = 0
    df_metadata = pd.read_csv('pairwise_contrasts_HVG1000_metadata_rep3.txt', sep = '\t', header = 0, index_col = 0)
    df_cluster = pd.read_csv('fda_clusters.txt', sep = '\t', header = 0, index_col = 0)

    fd, genes = build_fda(df_zscores, df_metadata, cell_type = 'cMono', perturbation = 'COVID_SEV')
    perform_smoothing(fd, df_cluster, genes, n_neighbors = 2, bandwidth = 1, 
                        cluster_key = 'clusters_fuzzy', output_folder = 'figures/smoothing_percluster/')