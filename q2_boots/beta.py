# ----------------------------------------------------------------------------
# Copyright (c) 2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from hdmedians import medoid
import pandas as pd
import numpy as np
import math


def beta_collection(ctx, table, metric, sampling_depth, phylogeny=None,
                    bypass_tips=False, n_threads=1, n=1000, random_seed=None,
                    with_replacement=True, pseudocount=1, alpha=None,
                    variance_adjusted=False):

    _resample = ctx.get_action('boots', 'resample')
    _beta = ctx.get_action('diversity', 'beta')
    _beta_phylogenetic = ctx.get_action('diversity', 'beta_phylogenetic')

    tables = _resample(table=table, sampling_depth=sampling_depth, n=n,
                       random_seed=random_seed, with_replacement=with_replacement)

    diversified_tables = []

    if phylogeny is not None:
        for passthrough in tables:
            diversified_tables.append(_beta(table=passthrough, metric=metric,
                                      pseudocount=pseudocount, n_jobs=n_threads))
    else:
        for passthrough in tables:
            diversified_tables.append(_beta_phylogenetic(
                table=passthrough,
                phylogeny=phylogeny,
                metric=metric,
                threads=n_threads,
                variance_adjusted=variance_adjusted,
                alpha=alpha,
                bypass_tips=bypass_tips
            ))

    return diversified_tables


def beta(ctx, table, metric, sampling_depth, representative, phylogeny=None,
         bypass_tips=False, n_threads=1, n=1000, random_seed=None,
         with_replacement=True, pseudocount=1, alpha=None, variance_adjusted=False):

    _beta = ctx.get_action('boots', 'beta_collection')

    matrices = _beta(table=table,
                     phylogeny=phylogeny,
                     metric=metric,
                     sampling_depth=sampling_depth,
                     n=n,
                     threads=n_threads,
                     variance_adjusted=variance_adjusted,
                     alpha=alpha,
                     bypass_tips=bypass_tips,
                     random_seed=random_seed)

    if representative == 'medoid':
        return get_medoid(matrices)
    elif representative == 'non-metric-mean':
        return per_cell_average(matrices, 'mean')
    elif representative == 'non-metric-median':
        return per_cell_average(matrices, 'median')


def per_cell_average(a, representation):
    w, h = a[0].shape

    mtx = []

    for col in range(w):
        c_row = []
        for row in range(h):
            cell = pd.Series([x[row][col] for x in a])

            if representation == 'median':
                c_row.append(cell.median())

            elif representation == 'mean':
                c_row.append(cell.mean())

        mtx.append(c_row.copy())

    return pd.DataFrame(mtx)


def get_medoid(a):

    arrays = np.array([x.values.ravel() for x in a]).T
    rep = medoid(arrays)

    side = int(math.sqrt(len(rep)))
    rep = np.reshape(np.array(rep, type(rep[0])), (side, side))

    return pd.DataFrame(rep)
