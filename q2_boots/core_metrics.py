# ----------------------------------------------------------------------------
# Copyright (c) 2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from hdmedians import medoid
from q2_boots.beta import per_cell_average
import pandas as pd


def core_metrics(ctx, table, sampling_depth, metric, metadata,
                 n_jobs, phylogeny=None, n=1, alpha_method='median',
                 beta_method='medoid', with_replacement=True,
                 random_seed=None):

    bootstrap = ctx.get_action('boots', 'resample')
    observed_features = ctx.get_action("diversity_lib", "observed_features")
    pielou_e = ctx.get_action('diversity_lib', 'pielou_evenness')
    shannon = ctx.get_action('diversity_lib', 'shannon_entropy')
    braycurtis = ctx.get_action('diversity_lib', 'bray_curtis')
    jaccard = ctx.get_action('diversity_lib', 'jaccard')
    pcoa = ctx.get_action('diversity', 'pcoa')
    emperor_plot = ctx.get_action('emperor', 'plot')

    results = []
    bootstrapped_tables = bootstrap(table=table,
                                    sampling_depth=sampling_depth,
                                    n=n, with_replacement=with_replacement,
                                    random_seed=random_seed)

    for m in (observed_features, shannon, pielou_e):
        results += alpha_representative(m, bootstrapped_tables, alpha_method)

    dms = []
    for m in (jaccard, braycurtis):
        beta_results = beta_representative(m, bootstrapped_tables, beta_method)
        results += beta_results
        dms += beta_results

    pcoas = []
    for dm in dms:
        pcoa_results = pcoa(distance_matrix=dm)
        results += pcoa_results
        pcoas += pcoa_results

    for pcoa in pcoas:
        results += emperor_plot(pcoa=pcoa, metadata=metadata)

    return tuple(results)


def beta_representative(func, tables, method):
    metric = []
    for table in tables:
        metric.append(func(table))

    if method == 'medoid':
        return medoid(metric)

    return per_cell_average(metric, representation=method)


def alpha_representative(func, tables, method):

    alpha = pd.DataFrame()
    for table in tables:
        alpha.append(func(table=table))

    if method == 'median':
        return alpha.median()
    elif method == 'mean':
        return alpha.mean()
