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
        print('failed')

    l_train = torch.clamp(p_means, max = 1e8)
    counts = Poisson(l_train).sample()
    return counts


def create_simulated_data(
    N = 100, # number of cells per group
    M = 60, # number of genes
    n_types = 2, # number of cell types
    n_perts = 3, # number of perturbation groups (including ctrl)
    n_times = 40, # number of timepoints
    n_timepatterns = 3, # number of time patterns
    pattern_coeff_gap = 4.0, # coefficient difference between time patterns
    missing_rate = 0.2, # ratio of missing time points
    base_mean_loc = 4, # base mean expression levels of genes
    base_mean_scale = 0.3, # base variance of expression levels of genes
    mean_type = 0, # mean of cell-type coefficients
    scale_type = 0.4, # variance of cell-type coefficients
    mean_pert = 0, # mean of perturbation coefficients
    scale_pert = 0.4, # variance of perturbation coefficient
    mean_type_pert = 0, # mean of cell-type & perturbation interaction coefficients
    scale_type_pert = 0.5, # variance of cell-type & perturbation interaction coefficients
    theta = 10.0, # inverse dispersion of Negative Binomial Distribution
    noise_level = 1, # random noise level of expression (not implemented)
    random_seed = 42,
    add_timecovariate = True,
    add_noise = True
):
    times = np.linspace(0, 1, n_times)
    times = random.sample(list(times), int(n_times * (1-missing_rate))) # apply time sparsity
    times.sort()
    print(times)

    # initializing time covariates
    if add_timecovariate:
        n_genes_pattern = int(M / n_timepatterns)
        # example 
        # time_covariates = np.array([-4, 0.0, 4.0])
        # time_intercepts = np.array([4, 1.0, 0.0])

        time_covariates = np.arange(np.floor(n_timepatterns / 2) * (-pattern_coeff_gap), 
                                    np.floor(n_timepatterns / 2) * (-pattern_coeff_gap) + n_timepatterns * pattern_coeff_gap, 
                                    pattern_coeff_gap)
        time_intercepts = []
        for i in time_covariates:
            if i < 0: time_intercepts.append(-i)
            elif i == 0: time_intercepts.append(1.)
            else: time_intercepts.append(0.)

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
                ht_linear = lambda t: (time_covar * t + time_inter)
                bgcp[i * n_genes_pattern : (i + 1) * n_genes_pattern, :] = ht_linear(time) * bgcp[i * n_genes_pattern : (i + 1) * n_genes_pattern, :]
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
    adata = create_simulated_data(
        N = 100, # number of cells per group
        M = 60,
        n_types = 2, # number of cell types
        n_perts = 3, # number of perturbation groups (including ctrl)
        n_times = 40,
        n_timepatterns = 3,
        pattern_coeff_gap = 5.0,
        missing_rate = 0.2,
        base_mean_loc = 4,
        base_mean_scale = 0.01,
        mean_type = 0,
        scale_type = 0.01,
        mean_pert = 0.3,
        scale_pert = 0.01,
        mean_type_pert = 2.0,
        scale_type_pert = 0.01,
        theta = 10.0,
        noise_level = 1,
        random_seed = 41,
        add_timecovariate = True,
        add_noise = True,
    )
    # adata.write('test_data.h5ad')
    
    random_seed_init = 41
    n_reps = 3

    # 1. simulation - noise level
    if not os.path.isdir('simulation_noise_level'):
        os.mkdir('simulation_noise_level')
    os.chdir('simulation_noise_level')

    for idx in range(n_reps):
        adata = create_simulated_data(
            N = 100, M = 60, n_types = 2, n_perts = 3, n_times = 20, n_timepatterns = 3,
            pattern_coeff_gap = 1.0, missing_rate = 0, base_mean_loc = 2, base_mean_scale = 0.1,
            mean_type = 0, scale_type = 0.2, mean_pert = 0.3, scale_pert = 0.2, mean_type_pert = 2.0,
            scale_type_pert = 0.2, noise_level = 0, theta = 10.0, random_seed = random_seed_init + idx, add_timecovariate = True, add_noise = False
        )
        adata.write('simulation_noise_0_rep' + str(idx) + '.h5ad')

        noise_levels = np.arange(1, 10)
        for noise_level in noise_levels:
            adata = create_simulated_data(
                N = 100, M = 60, n_types = 2, n_perts = 3, n_times = 20, n_timepatterns = 3,
                pattern_coeff_gap = 1.0, missing_rate = 0, base_mean_loc = 2, base_mean_scale = 0.1,
                mean_type = 0, scale_type = 0.2, mean_pert = 0.3, scale_pert = 0.2, mean_type_pert = 2.0,
                scale_type_pert = 0.2, noise_level = noise_level, theta = 10.0, random_seed = random_seed_init + idx, add_timecovariate = True, add_noise = True
            )
            adata.write('simulation_noise_' + str(noise_level) + '_rep' + str(idx) + '.h5ad')

    os.chdir('../')

    # 2. simulation - num of time points
    if not os.path.isdir('simulation_n_times'):
        os.mkdir('simulation_n_times')
    os.chdir('simulation_n_times')

    times = np.arange(5, 50, 5)
    for idx in range(n_reps):
        for time in times:
            adata = create_simulated_data(
                N = 100, M = 60, n_types = 2, n_perts = 3, n_times = time, n_timepatterns = 3,
                pattern_coeff_gap = 1.0, missing_rate = 0, base_mean_loc = 2, base_mean_scale = 0.1,
                mean_type = 0, scale_type = 0.2, mean_pert = 0.3, scale_pert = 0.2, mean_type_pert = 2.0,
                scale_type_pert = 0.2, noise_level = 0, theta = 10.0, random_seed = random_seed_init + idx, add_timecovariate = True, add_noise = False
            )
            adata.write('simulation_n_times_' + str(time) + '_rep' + str(idx) + '.h5ad')
    os.chdir('../')
    
    # 3. simulation - ratio of missing time points
    if not os.path.isdir('simulation_missing_rate_'):
        os.mkdir('simulation_missing_rate_')
    os.chdir('simulation_missing_rate_')

    missing_ratios = np.arange(0, 0.8, 0.1)
    for idx in range(n_reps):
        for missing_ratio in missing_ratios:
            adata = create_simulated_data(
                N = 100, M = 60, n_types = 2, n_perts = 3, n_times = 20, n_timepatterns = 3,
                pattern_coeff_gap = 1.0, missing_rate = missing_ratio, base_mean_loc = 2, base_mean_scale = 0.1,
                mean_type = 0, scale_type = 0.2, mean_pert = 0.3, scale_pert = 0.2, mean_type_pert = 2.0,
                scale_type_pert = 0.2, noise_level = 0, theta = 10.0, random_seed = random_seed_init + idx, add_timecovariate = True, add_noise = False
            )
            adata.write('simulation_missing_rate_' + str(missing_ratio) + '_rep' + str(idx) + '.h5ad')
    os.chdir('../')

    # 4. simulation - sequencing depth
    if not os.path.isdir('simulation_sequencing_depth_'):
        os.mkdir('simulation_sequencing_depth_')
    os.chdir('simulation_sequencing_depth_')

    base_means = np.arange(0.4, 4, 0.3)
    for idx in range(n_reps):
        for base_mean in base_means:
            try:
                print(base_mean)
                adata = create_simulated_data(
                    N = 100, M = 60, n_types = 2, n_perts = 3, n_times = 20, n_timepatterns = 3,
                    pattern_coeff_gap = 1.0, missing_rate = 0, base_mean_loc = base_mean, base_mean_scale = 0.01,
                    mean_type = 0, scale_type = 0.2, mean_pert = 0.3, scale_pert = 0.2, mean_type_pert = 2.0,
                    scale_type_pert = 0.2, noise_level = 0, theta = 10.0, random_seed = random_seed_init + idx + 2, add_timecovariate = True, add_noise = False
                )
                adata.write('simulation_sequencing_depth_' + str(np.round(base_mean, 2)) + '_rep' + str(idx) + '.h5ad')
            except:
                continue
    os.chdir('../')

    # 5. simulation - pattern gap
    if not os.path.isdir('simulation_pattern_gap_'):
        os.mkdir('simulation_pattern_gap_')
    os.chdir('simulation_pattern_gap_')

    pattern_gaps = np.arange(0.5, 4, 0.3)
    for idx in range(n_reps):
        for pattern_gap in pattern_gaps:
            adata = create_simulated_data(
                N = 100, M = 60, n_types = 2, n_perts = 3, n_times = 20, n_timepatterns = 3,
                pattern_coeff_gap = pattern_gap, missing_rate = 0, base_mean_loc = 2.0, base_mean_scale = 0.01,
                mean_type = 0, scale_type = 0.2, mean_pert = 0.3, scale_pert = 0.2, mean_type_pert = 2.0,
                scale_type_pert = 0.2, noise_level = 0, theta = 10.0, random_seed = random_seed_init + idx, add_timecovariate = True, add_noise = False
            )
            adata.write('simulation_pattern_gap_' + str(pattern_gap) + '_rep' + str(idx) + '.h5ad')
    os.chdir('../')