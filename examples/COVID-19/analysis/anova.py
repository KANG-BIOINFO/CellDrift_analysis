# use anova between covid-mild vs. other severe diseases, and use distance to extract what we need
import os
import pandas as pd
import numpy as np
import pickle
from dtw import dtw

import matplotlib as mpl
import matplotlib.pyplot as plt

import skfda
from skfda.inference.anova import oneway_anova
from skfda.representation import FDataGrid, FDataBasis

# global parameters
perts = ['COVID_CRIT','Sepsis','Flu','COVID_SEV','COVID_MILD','COVID_HCW_MILD']
groups = ['infection_severe', 'infection_severe', 'infection_severe', 'infection_severe', 'infection_mild', 'infection_mild']
group_dict = dict(zip(perts, groups))

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['axes.titlesize'] ='x-large'
mpl.rcParams['font.family'] = 'Arial'
mpl.rcParams['font.family'] = 'Arial'

#  remove duplicate time alignments
def make_unique_alignments(aligned_zvalues, aligned_timepoint):
    # get average z scores for multiple query -> one ref alignment
    df_align = pd.DataFrame({'aligned_timepoint': aligned_timepoint, 'aligned_zvalues': aligned_zvalues})
    df_avg = pd.DataFrame(df_align.groupby(by = ['aligned_timepoint']).mean()['aligned_zvalues'])

    updated_timepoint = list(df_avg.index.values)
    updated_zvalues = list(df_avg['aligned_zvalues'])
    return (updated_zvalues, updated_timepoint)


def create_fd(df_zscore, df_meta, gene, n_reps = 3, cell_type = 'cMono', 
                ref_perturbation = 'COVID_SEV', perturbations = perts, 
                step_pattern = 'symmetric2', alternative_pattern = 'asymmetric'):
    # filtering
    df_meta = df_meta.loc[df_meta['cell_type'] == cell_type, :].copy()
    df_meta = df_meta.loc[df_meta['perturbation'].isin(perturbations), :].copy()

    df_zscore = df_zscore.loc[[gene], :].copy()
    
    # order
    data_matrix = []
    timepoints = []
    fd_groups = []
    samples = []

    perturbations = list(np.setdiff1d(perturbations, ref_perturbation))
    perturbations = [ref_perturbation] + perturbations

    for perturbation in perturbations:
        for rep in range(n_reps):
            # preparing data matrix
            fd_groups.append(group_dict[perturbation])
            samples.append(perturbation + '_rep' + str(rep))

            rep_name = 'rep' + str(rep + 1)
            df_meta_sub = df_meta.loc[(df_meta['perturbation'] == perturbation) & (df_meta['rep'] == rep_name), :].copy()
            df_meta_sub = df_meta_sub.sort_values(['time_v2'], ascending = True)
            comps = list(df_meta_sub.index.values)
            df_zscore_sub = df_zscore.loc[:, comps].copy()
            timepoint = df_meta_sub['time_v2'].to_numpy()

            # create raw data matrix
            zvalues = df_zscore_sub.loc[gene, :].values

            # ref time
            if perturbation == ref_perturbation:
                ref_zvalues = zvalues
                ref_timepoint = timepoint
                # print(len(ref_timepoint))

                data_matrix.append(ref_zvalues)
                timepoints.append(ref_timepoint)

            # query time
            else:
                try:
                    # print(step_pattern)
                    alignment = dtw(zvalues, ref_zvalues, keep_internals = False, step_pattern = step_pattern, window_type = 'none')
                except:
                    # print(alternative_pattern)
                    alignment = dtw(zvalues, ref_zvalues, keep_internals = False, step_pattern = alternative_pattern, window_type = 'none')

                aligned_zvalues = np.array(zvalues)[alignment.index1]
                aligned_timepoint = np.array(ref_timepoint)[alignment.index2] # it contains multiple query -> one ref alignment

                # calcuate avg values for duplicate times
                updated_zvalues, updated_timepoint = make_unique_alignments(aligned_zvalues, aligned_timepoint)
                # print(len(updated_timepoint))
                data_matrix.append(updated_zvalues)
                timepoints.append(updated_timepoint)

    # create representation of functional data
    fd = skfda.FDataGrid(
        data_matrix = data_matrix,
        grid_points = timepoints[0] # since all timepoints are identical, we only need one here.
    )

    # visualization
    fig, ax = plt.subplots()
    fd.plot(axes = ax)
    fig.savefig('figures/aligned_curves_' + gene + '.pdf', bbox_inches = 'tight')

    return fd, fd_groups


def anova_test(fd, fd_groups):
    mild_idx = [i for i in range(len(fd_groups)) if fd_groups[i] == 'infection_mild']
    severe_idx = [i for i in range(len(fd_groups)) if fd_groups[i] == 'infection_severe']

    fd_mild = fd[mild_idx]
    fd_severe = fd[severe_idx]

    mild_mean = np.mean([np.mean(fd_mild.data_matrix[i]) for i in range(fd_mild.shape[0])])
    severe_mean = np.mean([np.mean(fd_severe.data_matrix[i]) for i in range(fd_severe.shape[0])])
    gap = severe_mean - mild_mean

    val, pval = oneway_anova(fd_mild, fd_severe)
    return val, pval, gap

def measure_distance():
    1


if __name__ == '__main__':
    reps = [1, 2, 3]
    df_meta = pd.DataFrame()
    df_zscore = pd.DataFrame()

    for rep in reps:
        zscore_file = 'pairwise_zscores_HVG1000_combined_rep' + str(rep) + '.txt'
        meta_file = 'pairwise_contrasts_HVG1000_metadata_rep' + str(rep) + '.txt'    

        df_meta_rep = pd.read_csv(meta_file, sep = '\t', header = 0, index_col = 0)
        df_meta_rep['rep'] = 'rep' + str(rep)
        df_meta_rep['index'] = [(i + '_' + j) for (i, j) in zip (df_meta_rep.index.values, df_meta_rep['rep'])]
        df_meta_rep['time_v2'] = [int(i.split('d')[1]) for i in df_meta_rep['time']]
        df_meta_rep.set_index('index', inplace = True)
        
        df_zscore_rep = pd.read_csv(zscore_file, sep = '\t', header = 0, index_col = 0)
        df_zscore_rep.columns = [i + '_rep' + str(rep) for i in df_zscore_rep.columns]
        df_zscore_rep.values[np.isnan(df_zscore_rep.values)] = 0

        df_meta = pd.concat([df_meta, df_meta_rep], axis = 0)
        df_zscore = pd.concat([df_zscore, df_zscore_rep], axis = 1)

    df_cluster = pd.read_csv('../fpca/selected_genes_20_EMClust_clusters.txt', sep = '\t', header = 0, index_col = 0)
    
    # cluster_id = 9
    clusters_shown = [2, 11]
    for cluster_id in clusters_shown:
        genes = df_cluster.loc[df_cluster['clusters'] == cluster_id, 'names'].copy()
        # genes = ['S100A8', 'S100A9', 'S100A12', 'RPS11', 'RPS12', 'RPS13', 'IFIT3', 'IFI16', 'IFI35', 'IFI44', 'IFI44L', 'IFIT1', 'IFIT2']

        df_results = pd.DataFrame(columns = ['Statistic', 'p-value', 'mean-gap (severe - mild)'])
        for gene in genes:
            print(gene)
            fd, fd_groups = create_fd(df_zscore, df_meta, gene, n_reps = 3, ref_perturbation = 'COVID_SEV', cell_type = 'cMono', perturbations = perts)
            val, pval, mean_gap = anova_test(fd, fd_groups)
            df_results.loc[gene, :] = [val, pval, mean_gap]

        df_results.to_csv('avova_EMcluster_' + str(cluster_id) + '.txt', sep = '\t')