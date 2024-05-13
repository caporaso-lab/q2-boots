# ----------------------------------------------------------------------------
# Copyright (c) 2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from q2_diversity._alpha import METRICS as ALPHA_METRICS
from q2_diversity._beta import METRICS as BETA_METRICS


def core_metrics(ctx, table, sampling_depth, metadata,
                 n_jobs=1, phylogeny=None, n=100, alpha_method='median',
                 beta_method='non-metric-median', with_replacement=True,
                 random_seed=None, alpha_metrics=["observed_features",
                                                  "shannon",
                                                  "pielou_e",
                                                  "faith_pd"],
                 beta_metrics=["jaccard", "braycurtis",
                               "weighted_unifrac", "unweighted_unifrac"]):

    bootstrap = ctx.get_action('boots', 'resample')
    pcoa = ctx.get_action('diversity', 'pcoa')
    emperor_plot = ctx.get_action('emperor', 'plot')
    alpha_average = ctx.get_action('boots', 'alpha_average')
    beta_average = ctx.get_action('boots', 'beta_average')

    bootstrapped_tables, = bootstrap(table=table,
                                     sampling_depth=sampling_depth,
                                     n=n, with_replacement=with_replacement,
                                     random_seed=random_seed)

    tables = bootstrapped_tables.values()

    alpha_metric_names = alpha_metrics

    alpha_metrics = {'PHYLO': {}, 'NONPHYLO': {}}

    for metric in alpha_metric_names:
        if metric in ALPHA_METRICS["PHYLO"]["IMPL"] or metric in\
           ALPHA_METRICS["PHYLO"]["UNIMPL"]:
            metric = ALPHA_METRICS['NAME_TRANSLATIONS'][metric]
            alpha_metrics['PHYLO'][metric] = ctx.get_action('diversity_lib',
                                                            metric)
        else:
            if metric in ALPHA_METRICS['NONPHYLO']['IMPL']:
                metric = ALPHA_METRICS['NAME_TRANSLATIONS'][metric]
                alpha_metrics['NONPHYLO'][metric] = ctx.get_action('diversity_lib',
                                                                   metric)
            else:
                alpha_metrics['NONPHYLO'][metric] = ctx.get_action('diversity_lib',
                                                                   'alpha_passthrough')

    beta_metric_names = beta_metrics

    beta_metrics = {'PHYLO': {}, 'NONPHYLO': {}}

    for metric in beta_metric_names:
        if metric in BETA_METRICS["PHYLO"]["IMPL"] or metric in\
           BETA_METRICS["PHYLO"]["UNIMPL"]:
            if metric in ('unweighted_unifrac', 'weighted_unifrac'):
                metric = BETA_METRICS['NAME_TRANSLATIONS'][metric]
                beta_metrics['PHYLO'][metric] = ctx.get_action('diversity_lib', metric)
            else:
                # handle unimplemented unifracs
                beta_metrics['PHYLO'][metric] =\
                    ctx.get_action('diversity_lib', 'beta_phylogenetic_passthrough')
        else:

            if metric in BETA_METRICS['NONPHYLO']['IMPL']:
                metric = BETA_METRICS['NAME_TRANSLATIONS'][metric]
                beta_metrics['NONPHYLO'][metric] = ctx.get_action('diversity_lib',
                                                                  metric)
            else:
                beta_metrics['NONPHYLO'][metric] = ctx.get_action('diversity_lib',
                                                                  'beta_passthrough')

    if phylogeny is None and (len(alpha_metrics['PHYLO']) > 0 or
                              len(beta_metrics['PHYLO']) > 0):
        print("WARNING: NO PHYLOGENY PROVIDED, PHYLOGENIC METRICS WILL NOT BE RUN")

    alphas = {}
    for m in (alpha_metrics['NONPHYLO'].values()):
        alpha = alpha_representative(m, tables, alpha_method)
        res, = alpha_average(data=alpha, average_method=alpha_method)
        alphas[m.__name__] = res
    if phylogeny is not None:
        for m in (alpha_metrics['PHYLO'].values()):
            alpha = alpha_representative(m, tables, alpha_method, phylogeny=phylogeny)
            res, = alpha_average(data=alpha, average_method=alpha_method)
            alphas[m.__name__] = res

    dms = {}
    for m in (beta_metrics['NONPHYLO'].values()):
        beta_results = beta_representative(m, tables, beta_method)
        res, = beta_average(data=beta_results, representative=beta_method)
        dms[m.__name__] = res
    if phylogeny is not None:
        for m in (beta_metrics['PHYLO'].values()):
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
