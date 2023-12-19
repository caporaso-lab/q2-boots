# ----------------------------------------------------------------------------
# Copyright (c) 2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

def core_metrics(ctx, table, sampling_depth, metric, metadata,
                 n_jobs, phylogeny=None, n=1, average_method='median'):

    represent = ctx.get_action("boots", "alpha_bootstrap_representative")
    observed_features = ctx.get_action("diversity_lib", "observed_features")
    pielou_e = ctx.get_action('diversity_lib', 'pielou_evenness')
    shannon = ctx.get_action('diversity_lib', 'shannon_entropy')
    braycurtis = ctx.get_action('diversity_lib', 'bray_curtis')
    jaccard = ctx.get_action('diversity_lib', 'jaccard')
    pcoa = ctx.get_action('diversity', 'pcoa')
    emperor_plot = ctx.get_action('emperor', 'plot')

    results = []
    bootstrapped_table, = represent(table=table, sampling_depth=sampling_depth,
                                    phylogeny=phylogeny, metric=metric, n=n,
                                    average_method=average_method)

    results.append(bootstrapped_table)

    for m in (observed_features, shannon, pielou_e):
        results += m(table=bootstrapped_table)

    dms = []
    for m in (jaccard, braycurtis):
        beta_results = m(table=bootstrapped_table, n_jobs=n_jobs)
        results += beta_results
        dms += beta_results

    pcoas = []
    for dm in dms:
        pcoa_results = pcoa(distance_matrix=dm)
        results += pcoa_results
        pcoas += pcoa_results

    for pcoa in pcoas:
        results += emperor_plot(pcoa=pcoa, metadata=metadata)

    return tuple(results)


def core_metrics_phylogenic():
    pass
