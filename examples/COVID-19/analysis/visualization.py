import numpy as np
import pandas as pd
import matplotlib.pylab as plt

import skfda
import skfda.misc.hat_matrix as hm
import skfda.preprocessing.smoothing.validation as val
from skfda.preprocessing.smoothing import KernelSmoother

def perform_smoothing(df_zscores, df_metadata, genes, n_neighbors, cell_type = 'cMono', 
                        perturbation = 'COVID_SEV', seed = 42,
                        fig_name = 'figures_cluster'):

    # preparation
    df_metadata['time_v2'] = [int(i.split('d')[1]) for i in df_metadata['time']]
    df_metadata = df_metadata.loc[(df_metadata['cell_type'] == cell_type) & (df_metadata['perturbation'] == perturbation), :].copy()
    df_metadata = df_metadata.sort_values(['time_v2'], ascending = True)
    timepoints = list(df_metadata['time_v2'])

    comps = list(df_metadata.index.values)
    df_zscores = df_zscores[comps]
    df_zscores = df_zscores.loc[genes, :]

    data_matrix = []
    for gene in df_zscores.index.values:
        gene_zvalues = df_zscores.loc[gene, :]
        data_matrix.append(gene_zvalues)

    # create representation of functional data
    fd = skfda.FDataGrid(
        data_matrix = data_matrix,
        grid_points = timepoints
    )

    # test visual
    fig, axes = plt.subplots(1,2,figsize = (8,3))
    fd.plot(axes = axes[0])
    fd.scatter(s = 0.5, axes = axes[1])
    fig.savefig('figures/plot_line_scatters_S100.png', bbox_inches = 'tight')

    # do smoothing
    scale_factor = (
        (fd.domain_range[0][1] - fd.domain_range[0][0]) / len(fd.grid_points[0])
    )

    bandwidth = n_neighbors * scale_factor

    # K-nearest neighbours kernel smoothing.
    knn = val.SmoothingParameterSearch(
        KernelSmoother(kernel_estimator=hm.KNeighborsHatMatrix()),
        n_neighbors,
        param_name='kernel_estimator__n_neighbors',
    )
    knn.fit(fd)
    knn_fd = knn.transform(fd)

    # Local linear regression kernel smoothing.
    llr = val.SmoothingParameterSearch(
        KernelSmoother(kernel_estimator=hm.LocalLinearRegressionHatMatrix()),
        bandwidth,
        param_name='kernel_estimator__bandwidth',
    )
    llr.fit(fd)
    llr_fd = llr.transform(fd)

    # Nadaraya-Watson kernel smoothing.
    nw = val.SmoothingParameterSearch(
        KernelSmoother(kernel_estimator=hm.NadarayaWatsonHatMatrix()),
        bandwidth,
        param_name='kernel_estimator__bandwidth',
    )
    nw.fit(fd)
    nw_fd = nw.transform(fd)

    # visualize fitted curves
    fig, axes = plt.subplots(1,2,figsize = (8,3))
    knn_fd.plot(axes = axes[0])
    knn_fd.scatter(s = 0.5, axes = axes[1])
    fig.savefig('figures/plot_line_scatters_knn_smoothing_S100.png', bbox_inches = 'tight')    

    fig, axes = plt.subplots(1,2,figsize = (8,3))
    llr_fd.plot(axes = axes[0])
    llr_fd.scatter(s = 0.5, axes = axes[1])
    fig.savefig('figures/plot_line_scatters_LR_smoothing_S100.png', bbox_inches = 'tight')    

    fig, axes = plt.subplots(1,2,figsize = (8,3))
    nw_fd.plot(axes = axes[0])
    nw_fd.scatter(s = 0.5, axes = axes[1])
    fig.savefig('figures/plot_line_scatters_NW_smoothing_S100.png', bbox_inches = 'tight')    

if __name__ == '__main__':
    # initialization
    df_zscores = pd.read_csv('../factor_analysis/pairwise_zscores_HVG1000_combined_rep3.txt', sep = '\t', header = 0, index_col = 0)
    df_metadata = pd.read_csv('../factor_analysis/pairwise_contrasts_HVG1000_metadata_rep3.txt', sep = '\t', header = 0, index_col = 0)
    df_features = pd.read_csv('../fpca/selected_overlap_features.txt', sep = '\t', header = 0, index_col = 0)
    features = list(df_features['names'])
    df_zscores.values[np.isnan(df_zscores.values)] = 0

    # run clustering
    n_neighbors = np.arange(1,3)
    
    features = ['S100A8', 'S100A9', 'S100A12']
    perform_smoothing(df_zscores, df_metadata, features, n_neighbors = n_neighbors,
                    cell_type = 'cMono', perturbation = 'COVID_SEV', seed = 42)