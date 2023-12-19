# ----------------------------------------------------------------------------
# Copyright (c) 2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from hdmedians import medoid


def beta(ctx, table, metric, sampling_depth, pseudocount=1, n_jobs=1, n=1):

    _bootstrap = ctx.get_action('boots', 'bootstrap')
    _beta = ctx.get_action('diversity', 'beta')

    tables = _bootstrap(table=table, sampling_depth=sampling_depth, n=n)

    diversified_tables = []

    for passthough in tables:
        diversified_tables.append(_beta(table=passthough, metric=metric,
                                        pseudocount=pseudocount, n_jobs=n_jobs))

    return diversified_tables


def beta_phylogenetic(ctx,
                      table,
                      phylogeny,
                      metric,
                      sampling_depth,
                      n=1,
                      threads=1,
                      variance_adjusted=False,
                      alpha=None,
                      bypass_tips=False):

    _bootstrap = ctx.get_action('boots', 'bootstrap')
    _beta_phylogenetic = ctx.get_action('diversity', 'beta_phylogenetic')

    tables = _bootstrap(table=table, sampling_depth=sampling_depth, n=n)

    diversified_tables = []

    for passthrough in tables:
        diversified_tables.append(_beta_phylogenetic(
            table=passthrough,
            phylogeny=phylogeny,
            metric=metric,
            threads=threads,
            variance_adjusted=variance_adjusted,
            alpha=alpha,
            bypass_tips=bypass_tips
        ))

    return diversified_tables


def beta_representative(ctx, table, metric, sampling_depth, pseudocount=1,
                        n_jobs=1, n=1, representative='medoid'):
    _beta = ctx.get_action('boots', 'beta')

    matrices = _beta(table=table,
                     metric=metric,
                     sampling_depth=sampling_depth,
                     pseudocount=pseudocount,
                     n_jobs=n_jobs,
                     n=n)

    if representative == 'medoid':
        return medoid(matrices)
    elif representative == 'non-metric-mean':
        return None
    elif representative == 'non-metric-median':
        return None


def beta_phylogenetic_representative(ctx,
                                     table,
                                     phylogeny,
                                     metric,
                                     sampling_depth,
                                     n=1,
                                     threads=1,
                                     variance_adjusted=False,
                                     alpha=None,
                                     bypass_tips=False,
                                     representative='medoid'):

    _beta_phylogenetic = ctx.get_action('boots', 'beta_phylogenetic')

    matrices = _beta_phylogenetic(table=table,
                                  phylogeny=phylogeny,
                                  metric=metric,
                                  sampling_depth=sampling_depth,
                                  n=n,
                                  threads=threads,
                                  variance_adjusted=variance_adjusted,
                                  alpha=alpha,
                                  bypass_tips=bypass_tips)

    if representative == 'medoid':
        return medoid(matrices)
    elif representative == 'non-metric-mean':
        return None
    elif representative == 'non-metric-median':
        return None
