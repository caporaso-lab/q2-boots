# ----------------------------------------------------------------------------
# Copyright (c) 2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from qiime2.plugin import (Plugin, Int, Range, Collection,
                           Citations, Str, Choices, Bool, Float)

from q2_types.feature_table import (
    FeatureTable, Frequency, RelativeFrequency, PresenceAbsence
)
from q2_types.sample_data import AlphaDiversity, SampleData

from q2_types.tree import Phylogeny, Rooted

from q2_diversity_lib.alpha import METRICS as alpha_metrics
from q2_diversity_lib.beta import METRICS as beta_metrics
from q2_types.distance_matrix import DistanceMatrix

import q2_boots

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
    'Phylogenetic tree containing tip identifiers that correspond to the feature '
    'identifiers in the table. This tree can contain tip ids that are not present '
    'in the table, but all feature ids must be present in this tree.'
)

random_seed_description = (
    'A seed to allow multiple runs to have the same outcome if the same seed is '
    'included'
)

Citations = Citations.load('citations.bib', package='q2_boots')
plugin = Plugin(
    name='boots',
    version=q2_boots.__version__,
    website='https://github.com/qiime2/q2-boots',
    package='q2_boots',
    short_description='placeholder',
    description='placeholder'
)

plugin.pipelines.register_function(
    function=q2_boots.resample,
    inputs={'table': FeatureTable[Frequency]},
    parameters={'sampling_depth': Int % Range(1, None),
                'n': Int % Range(1, None),
                'with_replacement': Bool,
                'random_seed': Int},
    outputs={'subsampled_tables': Collection[FeatureTable[Frequency]]},
    input_descriptions={'table': 'The table to be subsampled'},
    parameter_descriptions={
        'sampling_depth': ('The total frequency that each sample should be '
                           'subsampled to. Samples where the sum of frequencies '
                           'is less than the sampling depth will be not be '
                           'included in the resulting table.'),
        'n': 'The number of times to subsample the input table.',
        'with_replacement': '',
        'random_seed': random_seed_description
    },
    output_descriptions={
        'subsampled_tables': 'A collection of n tables normalized to the specified '
                             'sampling depth'
    },
    name='Bootstrap',
    description='This pipeline is a repeated subsampling of a specified input table. '
                'N tables are produced normalized so the sum of each sample\'s '
                'frequency is equal to the sampling depth.'
)

plugin.pipelines.register_function(
    function=q2_boots.alpha_collection,
    inputs={'table': FeatureTable[Frequency | RelativeFrequency | PresenceAbsence],
            'phylogeny': Phylogeny[Rooted]},
    parameters={'sampling_depth': Int % Range(1, None),
                'metric': Str % Choices(alpha_metrics['NONPHYLO']['IMPL'] |
                                        alpha_metrics['NONPHYLO']['UNIMPL'] |
                                        alpha_metrics['PHYLO']['IMPL'] |
                                        alpha_metrics['PHYLO']['UNIMPL']),
                'n': Int % Range(1, None),
                'random_seed': Int},
    outputs={'sample_data': Collection[SampleData[AlphaDiversity]]},
    input_descriptions={'table': 'The table to be diversified',
                        'phylogeny': phylogeny_description},
    parameter_descriptions={
        'sampling_depth': ('The total frequency that each sample should be '
                           'subsampled to. Samples where the sum of frequencies '
                           'is less than the sampling depth will be not be '
                           'included in the resulting table.'),
        'metric': 'The alpha diversity metric to be computed.',
        'n': 'The number of times to subsample the input table.',
        'random_seed': random_seed_description
    },
    output_descriptions={
        'sample_data': 'A collection of Alpha Divsersity Sample Data',
    },
    name='Alpha Bootstrap',
    description='Subsamples the input table multiple times and provides those in a ' +
                'collection of feature tables of the same type as an input.'
)

plugin.pipelines.register_function(
    function=q2_boots.alpha,
    inputs={'table': FeatureTable[Frequency | RelativeFrequency | PresenceAbsence],
            'phylogeny': Phylogeny[Rooted]},
    parameters={'sampling_depth': Int % Range(1, None),
                'metric': Str % Choices(alpha_metrics['NONPHYLO']['IMPL'] |
                                        alpha_metrics['NONPHYLO']['UNIMPL'] |
                                        alpha_metrics['PHYLO']['IMPL'] |
                                        alpha_metrics['PHYLO']['UNIMPL']),
                'n': Int % Range(1, None),
                'average_method': Str % Choices(['median' , 'mean' , 'mode']),
                'random_seed': Int},
    outputs={'sample_data': SampleData[AlphaDiversity]},
    input_descriptions={'table': 'The table to be diversified',
                        'phylogeny': phylogeny_description},
    parameter_descriptions={
        'sampling_depth': ('The total frequency that each sample should be '
                           'subsampled to. Samples where the sum of frequencies '
                           'is less than the sampling depth will be not be '
                           'included in the resulting table.'),
        'metric': 'The alpha diversity metric to be computed.',
        'n': 'The number of times to subsample the input table.',
        'random_seed': random_seed_description
    },
    output_descriptions={
        'sample_data': 'Vector containing per-sample alpha diversities.',
    },
    name='Alpha Bootstrap Representative',
    description=''
)

plugin.pipelines.register_function(
    function=q2_boots.beta,
    inputs={'table':
            FeatureTable[Frequency | RelativeFrequency | PresenceAbsence],
            'phylogeny': Phylogeny[Rooted]},
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
                'with_replacement': Bool,
                'variance_adjusted': Bool,
                'representative': Str % Choices(['non-metric-mean',
                                                 'non-metric-median',
                                                 'medoid']),
                'alpha': Float % Range(0, 1, inclusive_end=True)},
    outputs=[('distance_matrix', DistanceMatrix)],
    input_descriptions={
        'table': ('The feature table containing the samples over which beta '
                  'diversity should be computed.')
    },
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

plugin.methods.register_function(
    function=q2_boots.alpha_average,
    inputs={
        'data': Collection[SampleData[AlphaDiversity]]
    },
    parameters={
        'average_method': Str % Choices('mean', 'median'),
    },
    outputs={
        'alpha_diversity': SampleData[AlphaDiversity]
    },
    name='Alpha Average',
    description='Average Alpha Collection'
)

plugin.pipelines.register_function(
    function=q2_boots.beta_collection,
    inputs={'table':
            FeatureTable[Frequency | RelativeFrequency | PresenceAbsence],
            'phylogeny': Phylogeny[Rooted]},
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
                'with_replacement': Bool,
                'variance_adjusted': Bool,
                'alpha': Float % Range(0, 1, inclusive_end=True)},
    outputs={
        'distance_matrices': Collection[DistanceMatrix]
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
    output_descriptions={
        'distance_matrix': 'representative distance matrix',
    },
    parameter_descriptions={
        'representative': 'The method by which the data is represented.' +
                          'Medoid currently only works with small datasets and, ' +
                          'more importantly, small n values.'
    },
    name='Beta Average',
    description='Average of a Collection of Distance Matrices'
)

# plugin.pipelines.register_function(
#     function=q2_boots.core_metrics,
#     inputs={
#         'table': FeatureTable[Frequency | RelativeFrequency |
#                               PresenceAbsence],
#         'phylogeny': Phylogeny[Rooted],
#     },
#     parameters={
#         'metadata': Metadata,
#         'n_jobs': Int % Range(1, None),
#         'n': Int % Range(1, None),
#         'sampling_depth': Int % Range(1, None),
#         'alpha_method': Str % Choices('mean', 'median'),
#         'beta_method': Str % Choices('non-metric-mean',
#                                      'non-metric-median',
#                                      'medoid'),
#         'with_replacement': Bool,
#         'random_seed': Int
#     },
#     outputs=[
#         ('rarefied_table', Collection[FeatureTable[Frequency]]),
#         ('alpha_diversity', Collection[SampleData[AlphaDiversity]]),
#         ('distance_matrices', Collection[DistanceMatrix]),
#         ('pcoas', Collection[PCoAResults]),
#         ('visualizations', Collection[Visualization]),
#     ],
#     output_descriptions={
#
#     },
#     name='Core Metrics',
#     description='Bootstrapped Core Metrics'
# )
