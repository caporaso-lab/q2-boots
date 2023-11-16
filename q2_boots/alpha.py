# ----------------------------------------------------------------------------
# Copyright (c) 2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import pandas as pd

from q2_diversity_lib.alpha import METRICS


def alpha_bootstrap(ctx, table, sampling_depth, metric, phylogeny=None, n=1):

    if phylogeny is not None and metric in METRICS['NONPHYLO']:
        raise ValueError('You must use a phylogenic metric')

    elif phylogeny is None and metric in METRICS['PHYLO']:
        raise ValueError('You must use a non-phylogenic metric')

    _bootstrap = ctx.get_action("boots", "bootstrap")
    _alpha = ctx.get_action("divserity", "alpha")
    _alpha_phylogenetic = ctx.get_action("diversity", "alpha_phylogenetic")

    tables = _bootstrap(table=table, sampling_depth=sampling_depth, n=n)
    diversified_tables = []

    for table in tables:
        if phylogeny is not None:
            diversified_tables.append(_alpha_phylogenetic(
                table=table, metric=metric, phylogeny=phylogeny))
        else:
            diversified_tables.append(_alpha(
                table=table, metric=metric))

    return diversified_tables


def alpha_bootstrap_representative(ctx, table, sampling_depth, metric, phylogeny=None,
                                   n=1, average_method='median'):

    if phylogeny is not None and metric in METRICS['NONPHYLO']:
        raise ValueError('You must use a phylogenic metric when phylogeny is included.')

    elif phylogeny is None and metric in METRICS['PHYLO']:
        raise ValueError('You must use a non-phylogenic metric when no phylogeny is' +
                         'included.')

    _alpha_bootstrap = ctx.get_actions("boots", "alpha_bootstrap")
    sample_data = _alpha_bootstrap(table=table, sampling_depth=sampling_depth,
                                   phylogeny=phylogeny, metric=metric, n=n)

    representative_sample_data = pd.DataFrame(sample_data)

    if average_method == "median":
        representative_sample_data = representative_sample_data.median(axis=1)
    elif average_method == "mean":
        representative_sample_data = representative_sample_data.mean(axis=1)
    elif average_method == 'mode':
        representative_sample_data = representative_sample_data.mode(axis=1)

    return representative_sample_data
