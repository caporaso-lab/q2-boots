# ----------------------------------------------------------------------------
# Copyright (c) 2024, Caporaso Lab (https://cap-lab.bio).
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import pandas as pd

from q2_diversity_lib.alpha import METRICS


def alpha_collection(ctx, table, sampling_depth, metric, phylogeny=None, n=1000,
                     random_seed=None, replacement=False):

    if phylogeny is None and (metric in METRICS['PHYLO']['IMPL'] or
                              metric in METRICS['PHYLO']['UNIMPL']):
        raise ValueError('You must use a non-phylogenic metric when no '
                         'phylogeny is included.')

    if phylogeny is not None and (metric in METRICS['NONPHYLO']['IMPL'] or
                                  metric in METRICS['NONPHYLO']['UNIMPL']):
        phylogeny = None

    _bootstrap = ctx.get_action("boots", "resample")
    _alpha = ctx.get_action("diversity", "alpha")
    _alpha_phylogenetic = ctx.get_action("diversity", "alpha_phylogenetic")

    tables, = _bootstrap(table=table, sampling_depth=sampling_depth, n=n,
                         random_seed=random_seed, replacement=replacement)
    diversified_tables = []

    for table in tables.values():
        if phylogeny is not None and \
           metric not in METRICS['PHYLO']['IMPL'] and \
           metric not in METRICS['PHYLO']['UNIMPL']:
            tmp, = _alpha_phylogenetic(table=table, metric=metric,
                                       phylogeny=phylogeny)
            diversified_tables.append(tmp)
        else:
            tmp, = _alpha(table=table, metric=metric)
            diversified_tables.append(tmp)

    return (diversified_tables)


def alpha(ctx, table, sampling_depth, metric, phylogeny=None,
          n=1, average_method='median', random_seed=None, replacement=False):

    _alpha_bootstrap = ctx.get_action("boots", "alpha_collection")
    _alpha_average = ctx.get_action('boots', 'alpha_average')
    sample_data, = _alpha_bootstrap(table=table, sampling_depth=sampling_depth,
                                    phylogeny=phylogeny, metric=metric, n=n,
                                    random_seed=random_seed,
                                    replacement=replacement)

    result, = _alpha_average(sample_data, average_method)

    return result


def alpha_average(data: pd.Series, average_method: str) -> pd.Series:

    tbl = None
    i = 0
    metric = ''

    for a in data.values():
        if tbl is None:
            metric = a.name
            a.name = i
            tbl = pd.DataFrame(a)
        else:
            a.name = i
            tbl.join(a)
        i += 1

    if average_method == "median":
        representative_sample_data = tbl.median(axis=1)
    elif average_method == "mean":
        representative_sample_data = tbl.mean(axis=1)
    elif average_method == 'mode':
        representative_sample_data = tbl.mode(axis=1)
    representative_sample_data.name = metric
    return representative_sample_data
