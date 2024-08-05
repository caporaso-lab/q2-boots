# ----------------------------------------------------------------------------
# Copyright (c) 2024, Caporaso Lab (https://cap-lab.bio).
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import pandas as pd

import qiime2

from q2_boots._alpha import (_validate_alpha_metric, _get_alpha_metric_action,
                             _alpha_collection_from_tables)
from q2_boots._beta import (_validate_beta_metric, _get_beta_metric_action,
                            _beta_collection_from_tables)


def core_metrics(ctx, table, sampling_depth, metadata, n, replacement,
                 n_jobs=1, phylogeny=None, alpha_average_method='median',
                 beta_average_method='non-metric-median'):

    resample_action = ctx.get_action('boots', 'resample')
    alpha_average_action = ctx.get_action('boots', 'alpha_average')
    beta_average_action = ctx.get_action('boots', 'beta_average')
    pcoa_action = ctx.get_action('diversity', 'pcoa')
    emperor_plot_action = ctx.get_action('emperor', 'plot')

    alpha_metrics = ['pielou_e', 'observed_features', 'shannon']
    beta_metrics = ['braycurtis', 'jaccard']
    if phylogeny is not None:
        alpha_metrics.append('faith_pd')
        beta_metrics.extend(['unweighted_unifrac', 'weighted_unifrac'])
    # this validation step is unnecessary right now, but sets the stage for
    # user-provided metrics
    for alpha_metric in alpha_metrics:
        _validate_alpha_metric(alpha_metric, phylogeny)
    for beta_metric in beta_metrics:
        _validate_beta_metric(beta_metric, phylogeny)

    resampled_tables, = resample_action(table=table,
                                        sampling_depth=sampling_depth,
                                        n=n,
                                        replacement=replacement)

    alpha_vectors = {}
    for alpha_metric in alpha_metrics:
        alpha_metric_action = _get_alpha_metric_action(
            ctx, alpha_metric, phylogeny)
        alpha_collection = _alpha_collection_from_tables(
            resampled_tables, alpha_metric_action)
        avg_alpha_vector, = alpha_average_action(
            alpha_collection, alpha_average_method)
        alpha_vectors[alpha_metric] = avg_alpha_vector

    beta_dms = {}
    for beta_metric in beta_metrics:
        beta_metric_action = _get_beta_metric_action(
            ctx, beta_metric, phylogeny)
        beta_collection = _beta_collection_from_tables(
            resampled_tables, beta_metric_action)
        avg_beta_dm, = beta_average_action(
            beta_collection, beta_average_method)
        beta_dms[beta_metric] = avg_beta_dm

    pcoas = {}
    emperor_plots = {}
    for key, dm in beta_dms.items():
        pcoa_results, = pcoa_action(dm)
        pcoas[key] = pcoa_results
        emperor_plots[key] = emperor_plot_action(pcoa=pcoa_results,
                                                 metadata=metadata)[0]

    return resampled_tables, alpha_vectors, beta_dms, pcoas, emperor_plots
