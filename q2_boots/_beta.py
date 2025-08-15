# ----------------------------------------------------------------------------
# Copyright (c) 2024, Caporaso Lab (https://cap-lab.bio).
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import functools
import warnings

import numpy as np
import skbio
from hdmedians import medoid
import kmedoids

from q2_diversity_lib.beta import METRICS

_METRIC_MOD_DEFAULTS = {
    'bypass_tips': False,
    'pseudocount': 1,
    'alpha': None,
    'variance_adjusted': False
}


def beta_average(data: skbio.DistanceMatrix,
                 average_method: str) -> skbio.DistanceMatrix:
    # I need to be able to index into data.values(). Come up with a more
    # efficient way to do this.
    data = list(data.values())
    if average_method == 'medoid':
        warnings.warn('The current implementation of medoid may require '
                      'prohibitively large amounts of memory for large data '
                      'sets or large values of `n`. You can track progress '
                      'on this issue here: '
                      'https://github.com/caporaso-lab/q2-boots/issues/10')
        return _medoid(data)
    elif average_method == 'non-metric-mean':
        return _per_cell_average(data, 'mean')
    elif average_method == 'non-metric-median':
        return _per_cell_average(data, 'median')
    else:
        raise ValueError(f"Unknown average method {average_method}. "
                         "Available options are non-metric-median, non-metric-"
                         "mean, and medoid.")


def beta_collection(
        ctx, table, metric, sampling_depth, n, replacement,
        phylogeny=None,
        bypass_tips=_METRIC_MOD_DEFAULTS['bypass_tips'],
        pseudocount=_METRIC_MOD_DEFAULTS['pseudocount'],
        alpha=_METRIC_MOD_DEFAULTS['alpha'],
        variance_adjusted=_METRIC_MOD_DEFAULTS['variance_adjusted']):
    _validate_beta_metric(metric, phylogeny)

    resample_action = ctx.get_action("boots", "resample")
    beta_metric_action = _get_beta_metric_action(
        ctx, metric, phylogeny, bypass_tips, pseudocount, alpha,
        variance_adjusted)

    tables, = resample_action(table=table,
                              sampling_depth=sampling_depth,
                              n=n,
                              replacement=replacement)
    results = _beta_collection_from_tables(tables, beta_metric_action)

    return results


def beta(ctx, table, metric, sampling_depth, n, replacement,
         average_method='non-metric-median', phylogeny=None,
         bypass_tips=_METRIC_MOD_DEFAULTS['bypass_tips'],
         pseudocount=_METRIC_MOD_DEFAULTS['pseudocount'],
         alpha=_METRIC_MOD_DEFAULTS['alpha'],
         variance_adjusted=_METRIC_MOD_DEFAULTS['variance_adjusted']):
    beta_collection_action = ctx.get_action('boots', 'beta_collection')
    beta_average_action = ctx.get_action('boots', 'beta_average')
    dms, = beta_collection_action(table=table,
                                  phylogeny=phylogeny,
                                  metric=metric,
                                  sampling_depth=sampling_depth,
                                  n=n,
                                  pseudocount=pseudocount,
                                  replacement=replacement,
                                  variance_adjusted=variance_adjusted,
                                  alpha=alpha,
                                  bypass_tips=bypass_tips)

    result, = beta_average_action(dms, average_method)
    return result


def _per_cell_average(a, average_method):
    shape = a[0].shape
    ids = a[0].ids
    if average_method == 'median':
        average_fn = np.median
    elif average_method == 'mean':
        average_fn = np.mean
    else:
        raise ValueError(f"Unknown average method {average_method}. "
                         "Available options are median and mean.")
    condensed_dms = np.asarray([dm.condensed_form() for dm in a])
    average_condensed_dm = average_fn(condensed_dms, axis=0)

    result = np.zeros(shape)
    dm_indices = np.triu_indices(shape[0], 1)
    result[(dm_indices[0], dm_indices[1])] = average_condensed_dm
    result[(dm_indices[1], dm_indices[0])] = average_condensed_dm

    return skbio.DistanceMatrix(result, ids=ids)


def _medoid(a):
    # i dont think I need to condense for kmedoids
    condensed_dms = np.asarray([dm.condensed_form() for dm in a])
    medoid_dm_index = kmedoids.fasterpam(condensed_dms, medoids=1)
    print(medoid_dm_index)
    return a[medoid_dm_index.medoids[0]]


def _validate_beta_metric(metric, phylogeny):
    if _is_phylogenetic_beta_metric(metric) and phylogeny is None:
        raise ValueError(f'Metric {metric} requires a phylogenetic tree.')


def _get_beta_metric_action(
        ctx, metric, phylogeny,
        bypass_tips=_METRIC_MOD_DEFAULTS['bypass_tips'],
        pseudocount=_METRIC_MOD_DEFAULTS['pseudocount'],
        alpha=_METRIC_MOD_DEFAULTS['alpha'],
        variance_adjusted=_METRIC_MOD_DEFAULTS['variance_adjusted']):
    if _is_phylogenetic_beta_metric(metric):
        beta_metric_action = ctx.get_action("diversity", "beta_phylogenetic")
        beta_metric_action = functools.partial(
            beta_metric_action,
            phylogeny=phylogeny,
            metric=metric,
            bypass_tips=bypass_tips,
            alpha=alpha,
            variance_adjusted=variance_adjusted)
    else:
        beta_metric_action = ctx.get_action("diversity", "beta")
        beta_metric_action = functools.partial(beta_metric_action,
                                               metric=metric,
                                               pseudocount=pseudocount)
    return beta_metric_action


def _is_phylogenetic_beta_metric(metric):
    return metric in METRICS['PHYLO']['IMPL'] | METRICS['PHYLO']['UNIMPL']


def _beta_collection_from_tables(tables, beta_metric_action):
    results = []
    for table in tables.values():
        results.append(beta_metric_action(table=table)[0])
    return results
