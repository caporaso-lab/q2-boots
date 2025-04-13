# ----------------------------------------------------------------------------
# Copyright (c) 2024, Caporaso Lab (https://cap-lab.bio).
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from skbio import OrdinationResults
from qiime2 import Metadata
from q2_boots._alpha import (_validate_alpha_metric, _get_alpha_metric_action,
                             _alpha_collection_from_tables)
from q2_boots._beta import (_validate_beta_metric, _get_beta_metric_action,
                            _beta_collection_from_tables)


def kmer_diversity(ctx, table, sequences, sampling_depth, metadata, n,
                   replacement, kmer_size=16, tfidf=False, max_df=1.0, min_df=1,
                   max_features=None, alpha_average_method='median',
                   beta_average_method='non-metric-median', pc_dimensions=3,
                   color_by=None, norm='None',
                   alpha_metrics=['pielou_e', 'observed_features', 'shannon'],
                   beta_metrics=['braycurtis', 'jaccard']):

    resample_action = ctx.get_action('boots', 'resample')
    kmerize_action = ctx.get_action('kmerizer', 'seqs_to_kmers')
    alpha_average_action = ctx.get_action('boots', 'alpha_average')
    beta_average_action = ctx.get_action('boots', 'beta_average')
    pcoa_action = ctx.get_action('diversity', 'pcoa')
    scatter_action = ctx.get_action('vizard', 'scatterplot_2d')

    for alpha_metric in alpha_metrics:
        _validate_alpha_metric(alpha_metric, phylogeny=None)
    for beta_metric in beta_metrics:
        _validate_beta_metric(beta_metric, phylogeny=None)

    resampled_tables, = resample_action(table=table,
                                        sampling_depth=sampling_depth,
                                        n=n,
                                        replacement=replacement)
    kmer_tables = {}
    for key, resampled_table in resampled_tables.items():
        kmer_table, = kmerize_action(
            sequences, resampled_table, kmer_size, tfidf, max_df, min_df,
            max_features, norm)
        kmer_tables[key] = kmer_table

    alpha_vectors = {}
    for alpha_metric in alpha_metrics:
        alpha_metric_action = _get_alpha_metric_action(
            ctx, alpha_metric, phylogeny=None)
        alpha_collection = _alpha_collection_from_tables(
            kmer_tables, alpha_metric_action)
        avg_alpha_vector, = alpha_average_action(
            alpha_collection, alpha_average_method)
        alpha_vectors[alpha_metric] = avg_alpha_vector
        metadata = avg_alpha_vector.view(Metadata).merge(metadata)

    beta_dms = {}
    for beta_metric in beta_metrics:
        beta_metric_action = _get_beta_metric_action(
            ctx, beta_metric, phylogeny=None)
        beta_collection = _beta_collection_from_tables(
            kmer_tables, beta_metric_action)
        avg_beta_dm, = beta_average_action(
            beta_collection, beta_average_method)
        beta_dms[beta_metric] = avg_beta_dm

    pcoas = {}
    for key, dm in beta_dms.items():
        pcoa_results, = pcoa_action(dm)
        pcoas[key] = pcoa_results

    for pcoa, name in zip(pcoas.values(), beta_metrics):
        pc_result = pcoa.view(OrdinationResults)
        prop_explained = pc_result.proportion_explained[:pc_dimensions].values
        pc_result = pcoa.view(Metadata).to_dataframe().iloc[:, :pc_dimensions]
        pc_result.columns = ['{0} {1} ({2}%)'.format(name, c, int(p * 100)) for
                             c, p in zip(pc_result.columns, prop_explained)]
        metadata = Metadata(pc_result).merge(metadata)

    scatter_plot, = scatter_action(metadata=metadata, color_by=color_by)

    return (resampled_tables, kmer_tables, alpha_vectors,
            beta_dms, pcoas, scatter_plot)
