# ----------------------------------------------------------------------------
# Copyright (c) 2024, Caporaso Lab (https://cap-lab.bio).
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import skbio
import pandas as pd
import qiime2


def table_factory():
    table = pd.DataFrame(data=[[0, 1, 1],
                               [10, 10, 9],
                               [30, 20, 9],
                               [42, 42, 9]],
                         columns=['F1', 'F2', 'F3'],
                         index=['S1', 'S2', 'S3', 'S4'])
    return qiime2.Artifact.import_data(
        "FeatureTable[Frequency]", table, view_type=pd.DataFrame)


def sequences_factory():
    sequences = pd.Series(data=[skbio.DNA('GGACCCCTACGCCCATGGTAAACCGACTGGTCGTACGTGA'),  # noqa: E501
                                skbio.DNA('ACACGGACCTAAGAGCCGACCGCGTACAAAGGCGGGTACGTGCATTGGTTCCGGATCGCCCCGTACATCCGAAGAGCGTC'),  # noqa: E501
                                skbio.DNA('ACCCCGCCGGGTCATCATCATGCCAGCGACTACCA')],  # noqa: E501
                          index=['F1', 'F2', 'F3'])
    return qiime2.Artifact.import_data(
        "FeatureData[Sequence]", sequences, view_type=pd.Series)


def metadata_factory():
    metadata = pd.DataFrame(['odd', 'even', 'odd', 'even'],
                            index=['S1', 'S2', 'S3', 'S4'],
                            columns=['even-or-odd'])
    metadata.index.name = 'sample-id'
    metadata = qiime2.Metadata(metadata)
    return metadata


def _resample_bootstrap_example(use):
    table = use.init_artifact('table', table_factory)

    resampled_tables, = use.action(
        use.UsageAction(plugin_id='boots',
                        action_id='resample'),
        use.UsageInputs(table=table,
                        sampling_depth=20,
                        n=10,
                        replacement=True),
        use.UsageOutputNames(resampled_tables='bootstrapped_tables')
    )


def _resample_rarefaction_example(use):
    table = use.init_artifact('table', table_factory)

    resampled_tables, = use.action(
        use.UsageAction(plugin_id='boots',
                        action_id='resample'),
        use.UsageInputs(table=table,
                        sampling_depth=20,
                        n=10,
                        replacement=False),
        use.UsageOutputNames(resampled_tables='rarefaction_tables')
    )


def _alpha_bootstrap_example(use):
    table = use.init_artifact('table', table_factory)

    alpha_bootstrap, = use.action(
        use.UsageAction(plugin_id='boots',
                        action_id='alpha'),
        use.UsageInputs(table=table,
                        sampling_depth=20,
                        metric='observed_features',
                        n=10,
                        replacement=True,
                        average_method='median'),
        use.UsageOutputNames(
            average_alpha_diversity='observed_features_bootstrapped')
    )


def _alpha_rarefaction_example(use):
    table = use.init_artifact('table', table_factory)

    alpha_rarefaction, = use.action(
        use.UsageAction(plugin_id='boots',
                        action_id='alpha'),
        use.UsageInputs(table=table,
                        sampling_depth=20,
                        metric='observed_features',
                        n=10,
                        replacement=False,
                        average_method='median'),
        use.UsageOutputNames(
            average_alpha_diversity='observed_features_rarefaction')
    )


def _beta_bootstrap_example(use):
    table = use.init_artifact('table', table_factory)

    beta_bootstrap, = use.action(
        use.UsageAction(plugin_id='boots',
                        action_id='beta'),
        use.UsageInputs(table=table,
                        metric='braycurtis',
                        sampling_depth=20,
                        n=10,
                        replacement=True,
                        average_method='medoid'),
        use.UsageOutputNames(
            average_distance_matrix='braycurtis_bootstrapped')
    )


def _beta_rarefaction_example(use):
    table = use.init_artifact('table', table_factory)

    beta_rarefaction, = use.action(
        use.UsageAction(plugin_id='boots',
                        action_id='beta'),
        use.UsageInputs(table=table,
                        sampling_depth=20,
                        metric='braycurtis',
                        n=10,
                        replacement=False,
                        average_method='medoid'),
        use.UsageOutputNames(
            average_distance_matrix='braycurtis_rarefaction')
    )


def _core_metrics_bootstrap_example(use):
    table = use.init_artifact('table', table_factory)
    metadata = use.init_metadata('metadata', metadata_factory)

    # There seems to be a bug in the handling of this example, but I'm not yet
    # sure what it is. First, the Collection[Visualization] is not being handled
    # correctly. On the command line, the --o-emperor_plots value is
    # bootstrap_emperor_plots.qzv, though it is a directory and the expected
    # files are inside of it.
    # Additionally, setting the return value of use.action to `core_metrics, `
    # fails when calling `qiime dev refresh-cache`. As a result, the variable
    # I'm setting the return value to here differs from those used in the
    # previous examples. That, in turn, is triggering flake8 to complain about
    # unused variables. So, some stuff to unpack here 🥁, but the following
    # mostly works. (Queue *Hal fixing light bulb* video.)
    # This all goes for `_core_metrics_rarefaction_example` as well.
    core_metrics = use.action(  # noqa: F841
        use.UsageAction(plugin_id='boots',
                        action_id='core_metrics'),
        use.UsageInputs(table=table,
                        metadata=metadata,
                        sampling_depth=20,
                        n=10,
                        replacement=True,
                        alpha_average_method='median',
                        beta_average_method='medoid'),
        use.UsageOutputNames(
            resampled_tables='bootstrap_tables',
            alpha_diversities='bootstrap_alpha_diversities',
            distance_matrices='bootstrap_distance_matrices',
            pcoas='bootstrap_pcoas',
            emperor_plots='bootstrap_emperor_plots'
            )
    )


def _core_metrics_rarefaction_example(use):
    table = use.init_artifact('table', table_factory)
    metadata = use.init_metadata('metadata', metadata_factory)

    core_metrics = use.action(  # noqa: F841
        use.UsageAction(plugin_id='boots',
                        action_id='core_metrics'),
        use.UsageInputs(table=table,
                        metadata=metadata,
                        sampling_depth=20,
                        n=10,
                        replacement=False,
                        alpha_average_method='median',
                        beta_average_method='medoid'),
        use.UsageOutputNames(
            resampled_tables='rarefaction_tables',
            alpha_diversities='rarefaction_alpha_diversities',
            distance_matrices='rarefaction_distance_matrices',
            pcoas='rarefaction_pcoas',
            emperor_plots='rarefaction_emperor_plots'
            )
    )


def _kmer_diversity_bootstrap_example(use):
    table = use.init_artifact('table', table_factory)
    sequences = use.init_artifact('sequences', sequences_factory)
    metadata = use.init_metadata('metadata', metadata_factory)

    # There seems to be a bug in the handling of this example, but I'm not yet
    # sure what it is. First, the Collection[Visualization] is not being handled
    # correctly. On the command line, the --o-emperor_plots value is
    # bootstrap_emperor_plots.qzv, though it is a directory and the expected
    # files are inside of it.
    # Additionally, setting the return value of use.action to `core_metrics, `
    # fails when calling `qiime dev refresh-cache`. As a result, the variable
    # I'm setting the return value to here differs from those used in the
    # previous examples. That, in turn, is triggering flake8 to complain about
    # unused variables. So, some stuff to unpack here 🥁, but the following
    # mostly works. (Queue *Hal fixing light bulb* video.)
    # This all goes for `_core_metrics_rarefaction_example` as well.
    kmer_diversity = use.action(  # noqa: F841
        use.UsageAction(plugin_id='boots',
                        action_id='kmer_diversity'),
        use.UsageInputs(table=table,
                        sequences=sequences,
                        metadata=metadata,
                        sampling_depth=20,
                        n=10,
                        kmer_size=5,
                        replacement=True,
                        alpha_average_method='median',
                        beta_average_method='medoid'),
        use.UsageOutputNames(
            resampled_tables='bootstrap_tables',
            kmer_tables='kmer_tables',
            alpha_diversities='bootstrap_alpha_diversities',
            distance_matrices='bootstrap_distance_matrices',
            pcoas='bootstrap_pcoas',
            scatter_plot='scatter_plot'
            )
    )
