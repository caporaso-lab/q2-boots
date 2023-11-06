# ----------------------------------------------------------------------------
# Copyright (c) 2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------


def alpha_bootstrap(ctx, table, sampling_depth, metric, n=1):

    _bootstrap = ctx.get_action("boots", "bootstrap")
    _alpha = ctx.get_action("divserity_lib", "alpha_passthrough")

    tables = _bootstrap(ctx, table, sampling_depth, n)
    diversified_tables = []

    for table in tables:
        diversified_tables.append(_alpha(table, metric))

    return diversified_tables


def alpha_phylogenetic_bootstrap(ctx, table, sampling_depth, phylogeny, metric, n=1):

    _bootstrap = ctx.get_action("boots", "bootstrap")
    _alpha_phylogenetic = ctx.get_action("diversity", "alpha_phylogenetic")

    tables = _bootstrap(ctx, table, sampling_depth, n)
    diversified_tables = []

    for table in tables:
        diversified_tables.append(_alpha_phylogenetic(ctx, table, phylogeny, metric))

    return diversified_tables
