# ----------------------------------------------------------------------------
# Copyright (c) 2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import pandas as pd


def alpha_bootstrap(ctx, table, sampling_depth, metric, n=1):

    _bootstrap = ctx.get_action("boots", "bootstrap")
    _alpha = ctx.get_action("divserity", "alpha")

    tables = _bootstrap(table=table, sampling_depth=sampling_depth, n=n)
    diversified_tables = []

    for table in tables:
        diversified_tables.append(_alpha(table, metric))

    return diversified_tables


def alpha_phylogenetic_bootstrap(ctx, table, sampling_depth, phylogeny, metric, n=1):

    _bootstrap = ctx.get_action("boots", "bootstrap")
    _alpha_phylogenetic = ctx.get_action("diversity", "alpha_phylogenetic")

    tables = _bootstrap(table=table, sampling_depth=sampling_depth, n=n)
    diversified_tables = []

    for table in tables:
        diversified_tables.append(_alpha_phylogenetic(table=table, phylogeny=phylogeny,
                                                      metric=metric))

    return diversified_tables


def alpha_bootstrap_representative(ctx, table, sampling_depth, phylogeny, metric, n=1,
                                   average_method='median'):

    _alpha_bootstrap = ctx.get_actions("boots", "alpha_bootstrap")

    sample_data = _alpha_bootstrap(table=table, sampling_depth=sampling_depth,
                                   metric=metric, n=n)

    representative_sample_data = pd.DataFrame(sample_data)

    if average_method == "median":
        representative_sample_data = representative_sample_data.median(axis=1)
    elif average_method == "mean":
        representative_sample_data = representative_sample_data.mean(axis=1)

    return representative_sample_data
