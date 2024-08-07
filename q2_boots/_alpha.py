# ----------------------------------------------------------------------------
# Copyright (c) 2024, Caporaso Lab (https://cap-lab.bio).
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import functools

import pandas as pd

from q2_diversity import (alpha as q2div_alpha,
                          alpha_phylogenetic as q2div_alpha_phylogenetic)
from q2_diversity_lib.alpha import METRICS

from q2_boots._resample import resample


def alpha_average(data: pd.Series, average_method: str) -> pd.Series:
    if average_method == "median":
        average = pd.DataFrame.median
    elif average_method == "mean":
        average = pd.DataFrame.mean
    else:
        raise KeyError(f"Invalid average method: '{average_method}'. "
                       "Valid choices are 'median' and 'mean'.")
    data = pd.DataFrame(data.values())
    result = average(data, axis=0)
    result.name = data.index[0]
    return result


def alpha_collection(ctx, table, sampling_depth, metric, n,
                     replacement, phylogeny=None):
    _validate_alpha_metric(metric, phylogeny)

    resample_function = resample
    alpha_metric_function = _get_alpha_metric_function(ctx, metric, phylogeny)

    tables = resample_function(ctx=ctx,
                               table=table,
                               sampling_depth=sampling_depth,
                               n=n,
                               replacement=replacement)

    results = _alpha_collection_from_tables(tables, alpha_metric_function)

    return results


def alpha(ctx, table, sampling_depth, metric, n, replacement, phylogeny=None,
          average_method='median'):

    alpha_collection_function = alpha_collection
    alpha_average_action = ctx.get_action('boots', 'alpha_average')
    sample_data = alpha_collection_function(ctx=ctx,
                                            table=table,
                                            sampling_depth=sampling_depth,
                                            phylogeny=phylogeny,
                                            metric=metric,
                                            n=n,
                                            replacement=replacement)

    result, = alpha_average_action(sample_data, average_method)
    return result


def _validate_alpha_metric(metric, phylogeny):
    if _is_phylogenetic_alpha_metric(metric) and phylogeny is None:
        raise ValueError(f'Metric {metric} requires a phylogenetic tree.')


def _get_alpha_metric_function(ctx, metric, phylogeny):
    if _is_phylogenetic_alpha_metric(metric):
        alpha_metric_function = q2div_alpha_phylogenetic
        alpha_metric_function = functools.partial(alpha_metric_function,
                                                  ctx=ctx,
                                                  phylogeny=phylogeny,
                                                  metric=metric)
    else:
        alpha_metric_function = q2div_alpha
        alpha_metric_function = functools.partial(alpha_metric_function,
                                                  ctx=ctx,
                                                  metric=metric)
    return alpha_metric_function


def _is_phylogenetic_alpha_metric(metric):
    return metric in (METRICS['PHYLO']['IMPL'] | METRICS['PHYLO']['UNIMPL'])


def _alpha_collection_from_tables(tables, alpha_metric_function):
    results = []
    for table in tables.values():
        results.append(alpha_metric_function(table=table))
    return results
