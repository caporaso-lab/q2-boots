# ----------------------------------------------------------------------------
# Copyright (c) 2024, Caporaso Lab (https://cap-lab.bio).
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from qiime2.plugin import (Plugin, Int, Range, Collection,
                           Citations, Str, Choices, Bool, Float,
                           Metadata, Visualization)

from q2_types.feature_table import (
    FeatureTable, Frequency, RelativeFrequency, PresenceAbsence
)
from q2_types.sample_data import AlphaDiversity, SampleData

from q2_types.tree import Phylogeny, Rooted

from q2_diversity_lib.alpha import METRICS as alpha_metrics
from q2_diversity_lib.beta import METRICS as beta_metrics
from q2_types.distance_matrix import DistanceMatrix
from q2_types.ordination import PCoAResults

import q2_boots
from q2_boots._examples import (_resample_bootstrap_example,
                                _resample_rarefaction_example)

Citations = Citations.load('citations.bib', package='q2_boots')

plugin = Plugin(
    name='boots',
    version=q2_boots.__version__,
    website='https://github.com/caporaso-lab/q2-boots',
    package='q2_boots',
    short_description=('Bootstrapped and rarefaction-based diversity '
                       'analyses.'),
    description=('A plugin providing bootstrapped and rarefaction-based '
                 'diversity analyses, designed to mirror the interface of '
                 'q2-diversity.')
)


_feature_table_description = 'The feature table to be resampled.'
_sampling_depth_description = (
    'The total number of observations that each sample in `table` should be '
    'resampled to. Samples where the total number of observations in `table` '
    'is less than `sampling_depth` will be not be included in the output '
    'tables.')
_n_description = 'The number of resampled tables that should be generated.'
_replacement_description = (
    'Resample `table` with replacement (i.e., bootstrap) or without '
    'replacement (i.e., rarefaction).')
_resampled_tables_description = 'The `n` resampled tables.'

# Resampling

_resample_inputs = {
    'table': FeatureTable[Frequency]
}
_resample_parameters = {
    'sampling_depth': Int % Range(1, None),
    'n': Int % Range(1, None),
    'replacement': Bool
}
_resample_outputs = {
    'resampled_tables': Collection[FeatureTable[Frequency]]
}
_resample_input_descriptions = {
    'table': _feature_table_description
}
_resample_parameter_descriptions = {
    'sampling_depth': _sampling_depth_description,
    'n': _n_description,
    'replacement': _replacement_description
}
_resample_output_descriptions = {
    'resampled_tables': _resampled_tables_description
}

plugin.pipelines.register_function(
    function=q2_boots.resample,
    inputs=_resample_inputs,
    parameters=_resample_parameters,
    outputs=_resample_outputs,
    input_descriptions=_resample_input_descriptions,
    parameter_descriptions=_resample_parameter_descriptions,
    output_descriptions=_resample_output_descriptions,
    name='Resample feature table.',
    description=('Resample `table` to `sampling_depth` total observations with '
                 'replacement (i.e., bootstrapping) or without replacement '
                 '(i.e., rarefaction) `n` times, to generate `n` resampled '
                 'feature tables.'),
    examples={'Generate 10 bootstrapped tables.': _resample_bootstrap_example,
              'Generate 10 rarefied tables.': _resample_rarefaction_example}
)


n_jobs_description = (
    'The number of concurrent jobs to use in performing this calculation. '
    'May not exceed the number of available physical cores. If n_jobs = '
    '\'auto\', one job will be launched for each identified CPU core on the '
    'host.'
)

threads_description = (
    'The number of CPU threads to use in performing this calculation. '
    'May not exceed the number of available physical cores. If threads = '
    '\'auto\', one thread will be created for each identified CPU core on the '
    'host.'
)

phylogeny_description = (
    'Tree containing tip identifiers that correspond to the feature '
    'identifiers in the provided feature table. The tree can contain tips that '
    'are not present in the table, but all feature ids in the table must be '
    'present in this tree.'
)

random_seed_description = (
    'A seed to allow multiple runs to have the same outcome if the same seed '
    'is included.'
)

_alpha_inputs = {
    'table': FeatureTable[Frequency | RelativeFrequency | PresenceAbsence],
    'phylogeny': Phylogeny[Rooted]
}

_alpha_input_descriptions = {
    'table': 'The feature table to use in diversity computations.',
    'phylogeny': ('The phylogenetic tree to use in phylogenetic diversity '
                  'calculations. All feature ids in `table` must be present in '
                  'this tree, but this tree can contain feature ids that are '
                  'not present in `table`.')
}

_alpha_collection_parameters = {
    'sampling_depth': Int % Range(1, None),
    'metric': Str % Choices(alpha_metrics['NONPHYLO']['IMPL'] |
                            alpha_metrics['NONPHYLO']['UNIMPL'] |
                            alpha_metrics['PHYLO']['IMPL'] |
                            alpha_metrics['PHYLO']['UNIMPL']),
    'n': Int % Range(1, None),
    'replacement': Bool
}

_alpha_collection_parameter_descriptions = {
    'sampling_depth': _sampling_depth_description,
    'metric': 'The alpha diversity metric to be computed.',
    'n': _n_description,
    'replacement': _replacement_description
}

_alpha_average_parameters = {
    'average_method': Str % Choices('mean', 'median')
}

_alpha_average_parameter_descriptions = {
    'average_method': 'Method to use for averaging.'
}

plugin.methods.register_function(
    function=q2_boots.alpha_average,
    inputs={
        'data': Collection[SampleData[AlphaDiversity]]
    },
    parameters=_alpha_average_parameters,
    outputs={
        'averaged_alpha_diversities': SampleData[AlphaDiversity]
    },
    input_descriptions={
        'data': 'Alpha diversity vectors to be averaged.'
    },
    output_descriptions={
        'averaged_alpha_diversities':
        ('Per-sample alpha diversities representing the averages across the '
         'input data.')
    },
    parameter_descriptions=_alpha_average_parameter_descriptions,
    name='Average alpha diversity values.',
    description=('Compute the per-sample average across a collection of alpha '
                 'diversity vectors computed from the same samples.')
)

plugin.pipelines.register_function(
    function=q2_boots.alpha_collection,
    inputs=_alpha_inputs,
    parameters=_alpha_collection_parameters,
    outputs={'alpha_diversities': Collection[SampleData[AlphaDiversity]]},
    input_descriptions=_alpha_input_descriptions,
    parameter_descriptions=_alpha_collection_parameter_descriptions,
    output_descriptions={
        'alpha_diversities': ('`n` alpha diversity vectors, each containing '
                              'per-sample alpha diversity scores for the same '
                              'samples.'),
    },
    name='Iterative alpha diversity.',
    description=('Given a single feature table as input, this action resamples '
                 'the feature table `n` times to a total frequency of '
                 '`sampling depth` per sample, and then computes the specified '
                 'alpha diversity metric on each resulting `table`. The '
                 'resulting artifacts can be used, for example, to explore the '
                 'variance of `n` iterations of resampling.')
)

_alpha_parameters = _alpha_collection_parameters | _alpha_average_parameters
_alpha_parameter_descriptions = (_alpha_collection_parameter_descriptions |
                                 _alpha_average_parameter_descriptions)

plugin.pipelines.register_function(
    function=q2_boots.alpha,
    inputs=_alpha_inputs,
    parameters=_alpha_parameters,
    outputs={'sample_data': SampleData[AlphaDiversity]},
    input_descriptions=_alpha_input_descriptions,
    parameter_descriptions=_alpha_parameter_descriptions,
    output_descriptions={
        'sample_data': 'Vector containing per-sample alpha diversities.',
    },
    name='Bootstrapped or rarefaction-based alpha diversity.',
    description=('Given a single feature table as input, this action resamples '
                 'the feature table `n` times to a total frequency of '
                 '`sampling depth` per sample, and then computes the specified '
                 'alpha diversity metric on each resulting `table`. The '
                 'resulting artifacts are then averaged using the method '
                 'specified by `average_method`, and the resulting averaged '
                 'per-sample alpha diversities are returned.')
)

_beta_inputs = _alpha_inputs
_beta_input_descriptions = _alpha_input_descriptions

plugin.pipelines.register_function(
    function=q2_boots.beta,
    inputs=_beta_inputs,
    parameters={'metric': Str % Choices(beta_metrics['NONPHYLO']['IMPL'] |
                                        beta_metrics['NONPHYLO']['UNIMPL'] |
                                        beta_metrics['PHYLO']['IMPL'] |
                                        beta_metrics['PHYLO']['UNIMPL']),
                'pseudocount': Int % Range(1, None),
                'n_threads': Int % Range(1, None) | Str % Choices(['auto']),
                'n': Int % Range(1, None),
                'sampling_depth': Int % Range(1, None),
                'random_seed': Int,
                'bypass_tips': Bool,
                'replacement': Bool,
                'variance_adjusted': Bool,
                'representative': Str % Choices(['non-metric-mean',
                                                 'non-metric-median',
                                                 'medoid']),
                'alpha': Float % Range(0, 1, inclusive_end=True)},
    outputs=[('distance_matrix', DistanceMatrix)],
    input_descriptions=_beta_input_descriptions,
    parameter_descriptions={
        'metric': 'The beta diversity metric to be computed.',
        'pseudocount': ('A pseudocount to handle zeros for compositional '
                        'metrics.  This is ignored for other metrics.'),
        'random_seed': random_seed_description
    },
    output_descriptions={'distance_matrix': 'The resulting distance matrix.'},
    name='Beta diversity',
    description=("Computes a user-specified beta diversity metric for all "
                 "pairs of samples in a feature table.")
)

plugin.pipelines.register_function(
    function=q2_boots.beta_collection,
    inputs=_beta_inputs,
    parameters={'metric': Str % Choices(beta_metrics['NONPHYLO']['IMPL'] |
                                        beta_metrics['NONPHYLO']['UNIMPL'] |
                                        beta_metrics['PHYLO']['IMPL'] |
                                        beta_metrics['PHYLO']['UNIMPL']),
                'pseudocount': Int % Range(1, None),
                'replacement': Bool,
                'n_threads': Int % Range(1, None) | Str % Choices(['auto']),
                'n': Int % Range(1, None),
                'sampling_depth': Int % Range(1, None),
                'random_seed': Int,
                'bypass_tips': Bool,
                'replacement': Bool,
                'variance_adjusted': Bool,
                'alpha': Float % Range(0, 1, inclusive_end=True)},
    outputs={
        'distance_matrices': Collection[DistanceMatrix]
    },
    input_descriptions={

    },
    output_descriptions={

    },
    parameter_descriptions={

    },
    name='Beta Diversity Collection',
    description='Beta Diversity'
)

plugin.methods.register_function(
    function=q2_boots.beta_average,
    inputs={
        'data': Collection[DistanceMatrix],
    },
    parameters={
        'representative': Str % Choices(['non-metric-mean',
                                         'non-metric-median',
                                         'medoid']),
    },
    outputs=[
        ('distance_matrix', DistanceMatrix)
    ],
    input_descriptions={
        'data': 'Collection of Distance Matrices to be averaged'
    },
    output_descriptions={
        'distance_matrix': 'representative distance matrix',
    },
    parameter_descriptions={
        'representative': ''
    },
    name='Beta Average',
    description='Average of a Collection of Distance Matrices'
)

_core_metrics_inputs = _beta_inputs

plugin.pipelines.register_function(
    function=q2_boots.core_metrics,
    inputs=_core_metrics_inputs,
    parameters={
        'metadata': Metadata,
        'n_jobs': Int % Range(1, None),
        'n': Int % Range(1, None),
        'sampling_depth': Int % Range(1, None),
        'alpha_method': Str % Choices('mean', 'median'),
        'beta_method': Str % Choices('non-metric-mean',
                                     'non-metric-median',
                                     'medoid'),
        'replacement': Bool,
        'random_seed': Int
    },
    outputs=[
        ('rarefied_table', Collection[FeatureTable[Frequency]]),
        ('alpha_diversity', Collection[SampleData[AlphaDiversity]]),
        ('distance_matrices', Collection[DistanceMatrix]),
        ('pcoas', Collection[PCoAResults]),
        ('visualizations', Collection[Visualization]),
    ],
    input_descriptions={

    },
    parameter_descriptions={

    },
    output_descriptions={

    },
    name='Core Metrics',
    description='Bootstrapped Core Metrics'
)
