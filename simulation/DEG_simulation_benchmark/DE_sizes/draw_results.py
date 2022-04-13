import os
import pandas as pd 
import numpy as np 

import matplotlib as mpl 
import matplotlib.pyplot as plt 

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'Arial'

methods = ['CellDrift','t-test','wilcoxon','t-test-overdispersion','MAST']
colors = ['#EE3A3B', '#FACBCD', '#4D9B44', '#4177BC', '#B89BC9']
result_files = [i for i in os.listdir('./') if i.startswith('BenchDEfactor_metrics_threshold')]
for result_file in result_files:
    thred = result_file.split('_')[-1].split('.txt')[0]
    df_result = pd.read_csv(result_file, sep = '\t', header = 0, index_col = 0)

    de_factors = df_result['DEFacLoc'].unique()
    x_loc = np.arange(len(de_factors))
    
    for metric in ['TPR', 'FDR']:
        fig, ax = plt.subplots()

        for method, color in zip(methods,colors):
            df_sub = df_result.loc[df_result['Method'] == method, :].copy()  

            y_means = np.array(df_sub.groupby('DEFacLoc').mean()[metric])
            y_errs = np.array(df_sub.groupby('DEFacLoc').std()[metric])

            ax.plot(x_loc, y_means, color = color)
            ax.errorbar(x = x_loc, y = y_means, yerr = y_errs, color = color)
        fig.savefig('benchmark_results_thred_' + thred + '_' + metric + '.pdf', bbox_inches = 'tight')
            