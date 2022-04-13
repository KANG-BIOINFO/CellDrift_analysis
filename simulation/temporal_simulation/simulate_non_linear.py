import os
import torch
import anndata
import random
import numpy as np
import pandas as pd
from torch.distributions import Distribution, Gamma, Poisson, constraints
from torch.distributions.utils import broadcast_all

def _create_onehot_basic(N, c, p):
    rho_c = np.repeat( np.eye(c), N * p, axis = 0)
    rho_p = np.repeat( np.repeat( np.eye(p), N, axis = 0), c, axis = 0)
    rho_cp = np.repeat( np.eye(c * p), N, axis = 0)
    return rho_c, rho_p, rho_cp

def _gamma(theta, mu):
    concentration = theta
    rate = theta / mu
    # Important remark: Gamma is parametrized by the rate = 1/scale!
    gamma_d = Gamma(concentration=concentration, rate=rate)
    return gamma_d

def sample_NB(mu, theta):
    mu, theta = broadcast_all(mu, theta)
    try:
        gamma_d = _gamma(theta, mu)
        p_means = gamma_d.sample()
    except:
        print(torch.min(mu))
        print(torch.min(theta))
        print(mu)
        print(theta)

    l_train = torch.clamp(p_means, max = 1e8)
    counts = Poisson(l_train).sample()
    return counts


def create_simulated_data_nonlinear(
    N = 100, # number of cells per group
    M = 60, # number of genes
    n_types = 2, # number of cell types
    n_perts = 3, # number of perturbation groups (including ctrl)
    n_times = 20, # number of timepoints
    n_timepatterns = 3, # number of time patterns
    time_covariates = [-4.0, 0, 4.0], # time covariates
    time_intercepts = [1, 0.5, 0],
    time_center = 10,
    missing_rate = 0, # ratio of missing time points
    base_mean_loc = 2, # base mean expression levels of genes
    base_mean_scale = 0.1, # base variance of expression levels of genes
    mean_type = 0, # mean of cell-type coefficients
    scale_type = 0.2, # variance of cell-type coefficients
    mean_pert = 0.3, # mean of perturbation coefficients
    scale_pert = 0.2, # variance of perturbation coefficient
    mean_type_pert = 1.5, # mean of cell-type & perturbation interaction coefficients
    scale_type_pert = 0.2, # variance of cell-type & perturbation interaction coefficients
    theta = 10.0, # inverse dispersion of Negative Binomial Distribution
    noise_level = 1, # random noise level of expression (not implemented)
    random_seed = 42,
    add_timecovariate = True,
    add_noise = False
):
    times = np.linspace(0, 1, n_times)
    times = random.sample(list(times), int(n_times * (1-missing_rate))) # apply time sparsity
    times.sort()
    print(times)

    # initializing time covariates
    if add_timecovariate:
        n_genes_pattern = int(M / n_timepatterns)

        # example 
        # 0.1 * (x-10)^2 ; 5; -0.1 * (x-10)^2  + 10

    # generate count matrix for each timepoint
    for idx, time in enumerate(times):
        n_total = N * n_types * n_perts
        np.random.seed(random_seed)
        bg0 = np.random.normal(loc = base_mean_loc, scale = base_mean_scale, size = [1, M]) # intercepts for each gene without perturbation
        bg0 = np.repeat( bg0, n_total, axis = 0 )

        # covariates
        bgc = np.random.normal(loc = mean_type, scale = scale_type, size = [M, n_types])
        bgp = np.random.normal(loc = mean_pert, scale = scale_pert, size = [M, n_perts])
        bgp[:, 0] = 0
        bgcp = np.random.normal(loc = mean_type_pert, scale = scale_type_pert, size = [M, n_types * n_perts])
        bgcp[:, np.arange(0, n_types * n_perts, n_perts)] = 0 # make control groups as 0

        # add time covariates
        if add_timecovariate:
            for i in range(n_timepatterns):
                time_covar = time_covariates[i]
                time_inter = time_intercepts[i]
                # ht_linear = lambda t: (time_covar * t + time_inter) # linear function
                ht_nonlinear = lambda t: (time_covar * (time - time_center)**2 + time_inter)
                # bgcp[i * n_genes_pattern : (i + 1) * n_genes_pattern, :] = ht_linear(time) * bgcp[i * n_genes_pattern : (i + 1) * n_genes_pattern, :]
                bgcp[i * n_genes_pattern : (i + 1) * n_genes_pattern, :] = ht_nonlinear(time) * bgcp[i * n_genes_pattern : (i + 1) * n_genes_pattern, :]
        if add_noise:
            np.random.seed(random_seed + idx)
            noise = np.random.uniform(low = 0.0, high = 1.0, size = [M, n_types * n_perts]) / 10 * noise_level
            bgcp = bgcp + noise
            bgcp[:, np.arange(0, n_types * n_perts, n_perts)] = 0 # make control groups as 0

        # one-hot matrix
        rho_nc, rho_np, rho_ncp = _create_onehot_basic(N = N, c = n_types, p = n_perts)

        # shape should be n_cells x n_features
        mu = bg0 + np.dot(rho_nc, bgc.T) + np.dot(rho_np, bgp.T) + np.dot(rho_ncp, bgcp.T)

        # generate raw count matrix
        mu = torch.tensor(mu)
        theta = torch.tensor(theta)
        # print(mu[:5, :5])
        # print(theta)
        counts = sample_NB(mu = mu, theta = theta).numpy()

        ### transfer it to anndata
        barcodes = ['Cell_' + str(i) for i in range(n_total)]
        types = ['Type_' + str(i) for i in range(n_types)]
        types = np.repeat(types, N * n_perts, axis = 0)

        # set up obs
        perts = ['Control'] + ['Perturb_' + str(i) for i in range(n_perts - 1)]
        perts = list(np.repeat( perts, N, axis = 0 )) * n_types

        df_cells = pd.DataFrame({'barcodes': barcodes, 'cell_type': types, 'perturb': perts})
        df_cells = df_cells.set_index('barcodes')
        df_cells['time'] = time

        # set up var
        features = ['Gene_' + str(i) for i in range(M)]
        df_genes = pd.DataFrame(index = features)

        # add beta metadata
        types = ['Type_' + str(i) for i in range(n_types)]
        perts = ['Control'] + ['Perturb_' + str(i) for i in range(n_perts - 1)]
        types_perts = []
        for i in types:
            for j in perts:
                types_perts.append(i + '_' + j)
        var_cols = ['beta_' + i for i in types] + ['beta_' + i for i in perts] + ['beta_' + i for i in types_perts]
        var_cols = [(i + '_' + str(np.round(time, 2))) for i in var_cols]
        varm_sub = np.concatenate([bgc, bgp, bgcp], axis = 1)

        if idx == 0:
            varm = varm_sub.copy()
            beta_meta = var_cols
        else:
            varm = np.concatenate([varm, varm_sub], axis = 1)
            beta_meta += var_cols

        # set up anndata
        adata = anndata.AnnData(
            X = counts,
            obs = df_cells,
            var = df_genes
        )

        # concatenate
        if idx == 0:
            adata_combined = adata.copy()
        else:
            adata_combined = adata_combined.concatenate(adata, join = 'outer', index_unique = None)

    # write anndata
    adata_combined.varm['beta'] = varm

    adata_combined.uns['beta_metadata'] = beta_meta
    adata_combined.uns['theta'] = theta.numpy()
    adata_combined.uns['n_types'] = n_types
    adata_combined.uns['n_perts'] = n_perts
    adata_combined.uns['N'] = N
    adata_combined.uns['M'] = M

    adata_combined.obs['Cell'] = ['Cell_' + str(i) for i in range(adata_combined.shape[0])]
    adata_combined.obs.set_index('Cell', inplace = True)

    return adata_combined


if __name__ == '__main__':

    # adata = create_simulated_data(N = 100, M = 50, n_types = 2, n_perts = 3, n_times = 30, n_timepatterns = 3, random_seed = 41)
    time_covariates_list = [[-8.0, 0, 8.0], [-4.0, 0, 4.0], [-2.0, 0, 2.0]]
    time_intercepts_list = [[2.0, 1, 0], [1, 0.5, 0], [0.5, 0.25, 0]]
    covariate_levels = [1, 2, 3]
    n_reps = 3

    os.chdir('simulation_nonlinear_covariates')
    for rep_id in range(n_reps):
        for covariate_level, time_covariates, time_intercepts in zip(covariate_levels, time_covariates_list, time_intercepts_list):
            adata = create_simulated_data_nonlinear(
                N = 100, # number of cells per group
                M = 60, # number of genes
                n_types = 2, # number of cell types
                n_perts = 3, # number of perturbation groups (including ctrl)
                n_times = 20, # number of timepoints
                n_timepatterns = 3, # number of time patterns
                time_covariates = time_covariates, # time covariates
                time_intercepts = time_intercepts,
                time_center = 0.5,
                missing_rate = 0, # ratio of missing time points
                base_mean_loc = 2, # base mean expression levels of genes
                base_mean_scale = 0.1, # base variance of expression levels of genes
                mean_type = 0, # mean of cell-type coefficients
                scale_type = 0.2, # variance of cell-type coefficients
                mean_pert = 0.3, # mean of perturbation coefficients
                scale_pert = 0.2, # variance of perturbation coefficient
                mean_type_pert = 1.5, # mean of cell-type & perturbation interaction coefficients
                scale_type_pert = 0.2, # variance of cell-type & perturbation interaction coefficients
                theta = 10.0, # inverse dispersion of Negative Binomial Distribution
                noise_level = 1, # random noise level of expression (not implemented)
                random_seed = 42,
                add_timecovariate = True,
                add_noise = False
            )
            adata.write('simulation_nonlinear_covariate_level_' + str(covariate_level) + '_rep' + str(rep_id) +'.h5ad')