# ----------------------------------------------------------------------------
# Copyright (c) 2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

def core_metrics(ctx, table, sampling_depth, metadata,
                 n_jobs, phylogeny=None, n=1, alpha_method='median',
                 beta_method='non-metric-median', with_replacement=True,
                 random_seed=None):

    bootstrap = ctx.get_action('boots', 'resample')
    observed_features = ctx.get_action("diversity_lib", "observed_features")
    pielou_e = ctx.get_action('diversity_lib', 'pielou_evenness')
    shannon = ctx.get_action('diversity_lib', 'shannon_entropy')
    braycurtis = ctx.get_action('diversity_lib', 'bray_curtis')
    jaccard = ctx.get_action('diversity_lib', 'jaccard')
    pcoa = ctx.get_action('diversity', 'pcoa')
    emperor_plot = ctx.get_action('emperor', 'plot')
    alpha_average = ctx.get_action('boots', 'alpha_average')
    beta_average = ctx.get_action('boots', 'beta_average')
    faith_pd = ctx.get_action('diversity_lib', 'faith_pd')
    unweighted_unifrac = ctx.get_action('diversity_lib', 'unweighted_unifrac')
    weighted_unifrac = ctx.get_action('diversity_lib', 'weighted_unifrac')

    bootstrapped_tables, = bootstrap(table=table,
                                     sampling_depth=sampling_depth,
                                     n=n, with_replacement=with_replacement,
                                     random_seed=random_seed)

    tables = bootstrapped_tables.values()

    alphas = {}
    for m in (observed_features, shannon, pielou_e):
        alpha = alpha_representative(m, tables, alpha_method)
        res, = alpha_average(data=alpha, average_method=alpha_method)
        alphas[m.__name__] = res
    if phylogeny is not None:
        for m in (faith_pd, ):
            alpha = alpha_representative(m, tables, alpha_method, phylogeny=phylogeny)
            res, = alpha_average(data=alpha, average_method=alpha_method)
            alphas[m.__name__] = res

    dms = {}
    for m in (jaccard, braycurtis):
        beta_results = beta_representative(m, tables, beta_method)
        res, = beta_average(data=beta_results, representative=beta_method)
        dms[m.__name__] = res
    if phylogeny is not None:
        for m in (unweighted_unifrac, weighted_unifrac):
            beta_results = beta_representative(m, tables, beta_method, phylogeny,
                                               n_threads=n_jobs)
            beta_results, = beta_average(data=beta_results,
                                         representative=beta_method)
            dms[m.__name__] = beta_results

    pcoas = {}
    for key, dm in dms.items():
        pcoa_results, = pcoa(distance_matrix=dm)
        pcoas[key] = pcoa_results

    visualizations = {}
    for key, pcoa in pcoas.items():
        visualizations[key] = (emperor_plot(pcoa=pcoa, metadata=metadata)[0])

    return bootstrapped_tables, alphas, dms, pcoas, visualizations


def beta_representative(func, tables, method,
                        phylogeny=None, n_threads=1):
    metric = []
    if phylogeny is None:
        for table in tables:
            metric.append(func(table=table)[0])
    else:
        for table in tables:
            metric.append(func(table=table, phylogeny=phylogeny,
                               threads=n_threads)[0])

    return metric


def alpha_representative(func, tables, method, phylogeny=None):

    alpha = []
    if phylogeny is None:
        for table in tables:
            alpha.append(func(table=table)[0])
    else:
        for table in tables:
            alpha.append(func(table=table, phylogeny=phylogeny)[0])

    return alpha
