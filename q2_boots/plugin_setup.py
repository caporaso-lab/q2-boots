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

Citations = Citations.load('citations.bib', package='q2_boots')
plugin = Plugin(
    name='boots',
    version=q2_boots.__version__,
    website='https://github.com/qiime2/q2-boots',
    package='q2_boots',
    short_description='placeholder',
    description='placeholder'
)

plugin.methods.register_function(
    function=q2_boots._bootstrap_iteration,
    inputs={'table': FeatureTable[Frequency | RelativeFrequency |
                                  PresenceAbsence]},
    parameters={'sampling_depth': Int % Range(1, None)},
    outputs={'subsampled_table': FeatureTable[Frequency | RelativeFrequency |
                                              PresenceAbsence]},
    input_descriptions={'table': 'The table to be subsampled'},
    parameter_descriptions={
        'sampling_depth': ('The total frequency that each sample should be '
                           'subsampled to. Samples where the sum of frequencies '
                           'is less than the sampling depth will be not be '
                           'included in the resulting table.')},
    output_descriptions={'subsampled_table': 'The table that reaches the threshold '
                         'specified by sampling depth'},
    name='Bootstrap Iteration',
    description='This private function is the logic for each iteration of '
                'bootstrap'
)

plugin.pipelines.register_function(
    function=q2_boots.bootstrap,
    inputs={'table': FeatureTable[Frequency]},
    parameters={'sampling_depth': Int % Range(1, None),
                'n': Int % Range(1, None)},
    outputs={'subsampled_tables': Collection[FeatureTable[Frequency |
                                                          RelativeFrequency |
                                                          PresenceAbsence]]},
    input_descriptions={'table': 'The table to be subsampled'},
    parameter_descriptions={
        'sampling_depth': ('The total frequency that each sample should be '
                           'subsampled to. Samples where the sum of frequencies '
                           'is less than the sampling depth will be not be '
                           'included in the resulting table.'),
        'n': 'The number of times to subsample the input table.'
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
    function=q2_boots.alpha,
    inputs={'table': FeatureTable[Frequency | RelativeFrequency | PresenceAbsence],
            'phylogeny': Phylogeny[Rooted]},
    parameters={'sampling_depth': Int % Range(1, None),
                'metric': Str % Choices(alpha_metrics['NONPHYLO']['IMPL'] |
                                        alpha_metrics['NONPHYLO']['UNIMPL'] |
                                        alpha_metrics['PHYLO']['IMPL'] |
                                        alpha_metrics['PHYLO']['UNIMPL']),
                'n': Int % Range(1, None)},
    outputs={'sample_data': Collection[SampleData[AlphaDiversity]]},
    input_descriptions={'table': 'The table to be diversified',
                        'phylogeny': phylogeny_description},
    parameter_descriptions={
        'sampling_depth': ('The total frequency that each sample should be '
                           'subsampled to. Samples where the sum of frequencies '
                           'is less than the sampling depth will be not be '
                           'included in the resulting table.'),
        'metric': '',
        'n': 'The number of times to subsample the input table.'
    },
    output_descriptions={
        'sample_data': '',
    },
    name='Alpha Bootstrap',
    description='Subsamples the input table multiple times and provides those in a ' +
                'collection of feature tables of the same type as an input.'
)

plugin.pipelines.register_function(
    function=q2_boots.alpha_representative,
    inputs={'table': FeatureTable[Frequency | RelativeFrequency | PresenceAbsence],
            'phylogeny': Phylogeny[Rooted]},
    parameters={'sampling_depth': Int % Range(1, None),
                'metric': Str % Choices(alpha_metrics['NONPHYLO']['IMPL'] |
                                        alpha_metrics['NONPHYLO']['UNIMPL'] |
                                        alpha_metrics['PHYLO']['IMPL'] |
                                        alpha_metrics['PHYLO']['UNIMPL']),
                'n': Int % Range(1, None),
                'average_method': Str % Choices(['median' , 'mean' , 'mode'])},
    outputs={'sample_data': SampleData[AlphaDiversity]},
    input_descriptions={'table': 'The table to be diversified',
                        'phylogeny': phylogeny_description},
    parameter_descriptions={
        'sampling_depth': ('The total frequency that each sample should be '
                           'subsampled to. Samples where the sum of frequencies '
                           'is less than the sampling depth will be not be '
                           'included in the resulting table.'),
        'metric': '',
        'n': 'The number of times to subsample the input table.'
    },
    output_descriptions={
        'sample_data': '',
    },
    name='Alpha Bootstrap Representative',
    description=''
)

plugin.pipelines.register_function(
    function=q2_boots.beta,
    inputs={'table':
            FeatureTable[Frequency | RelativeFrequency | PresenceAbsence]},
    parameters={'metric': Str % Choices(beta_metrics['NONPHYLO']['IMPL'] |
                                        beta_metrics['NONPHYLO']['UNIMPL']),
                'pseudocount': Int % Range(1, None),
                'n_jobs': Int % Range(1, None) | Str % Choices(['auto']),
                'n': Int % Range(1, None),
                'sampling_depth': Int % Range(1, None)},
    outputs=[('distance_matrix', DistanceMatrix)],
    input_descriptions={
        'table': ('The feature table containing the samples over which beta '
                  'diversity should be computed.')
    },
    parameter_descriptions={
        'metric': 'The beta diversity metric to be computed.',
        'pseudocount': ('A pseudocount to handle zeros for compositional '
                        'metrics.  This is ignored for other metrics.'),
        'n_jobs': n_jobs_description
    },
    output_descriptions={'distance_matrix': 'The resulting distance matrix.'},
    name='Beta diversity',
    description=("Computes a user-specified beta diversity metric for all "
                 "pairs of samples in a feature table.")
)

plugin.pipelines.register_function(
    function=q2_boots.beta_phylogenetic,
    inputs={'table':
            FeatureTable[Frequency | RelativeFrequency | PresenceAbsence],
            'phylogeny': Phylogeny[Rooted]},
    parameters={'metric': Str % Choices(beta_metrics['PHYLO']['IMPL'] |
                                        beta_metrics['PHYLO']['UNIMPL']),
                'threads': Int % Range(1, None) | Str % Choices(['auto']),
                'variance_adjusted': Bool,
                'alpha': Float % Range(0, 1, inclusive_end=True),
                'bypass_tips': Bool,
                'n': Int % Range(1, None),
                'sampling_depth': Int % Range(1, None)},
    outputs=[('distance_matrix', DistanceMatrix)],
    input_descriptions={
        'table': ('The feature table containing the samples over which beta '
                  'diversity should be computed.'),
        'phylogeny': phylogeny_description,
    },
    parameter_descriptions={
        'metric': 'The beta diversity metric to be computed.',
        'threads': threads_description,
        'variance_adjusted': ('Perform variance adjustment based on Chang et '
                              'al. BMC Bioinformatics 2011. Weights distances '
                              'based on the proportion of the relative '
                              'abundance represented between the samples at a'
                              ' given node under evaluation.'),
        'alpha': ('This parameter is only used when the choice of metric is '
                  'generalized_unifrac. The value of alpha controls importance'
                  ' of sample proportions. 1.0 is weighted normalized UniFrac.'
                  ' 0.0 is close to unweighted UniFrac, but only if the sample'
                  ' proportions are dichotomized.'),
        'bypass_tips': ('In a bifurcating tree, the tips make up about 50% of '
                        'the nodes in a tree. By ignoring them, specificity '
                        'can be traded for reduced compute time. This has the'
                        ' effect of collapsing the phylogeny, and is analogous'
                        ' (in concept) to moving from 99% to 97% OTUs')
    },
    output_descriptions={'distance_matrix': 'The resulting distance matrix.'},
    name='Beta diversity (phylogenetic)',
    description=("Computes a user-specified phylogenetic beta diversity metric"
                 " for all pairs of samples in a feature table.")
)

plugin.pipelines.register_function(
    function=q2_boots.core_metrics,
    inputs={'table':
            FeatureTable[Frequency | RelativeFrequency | PresenceAbsence],
            'phylogeny': Phylogeny[Rooted]},
    parameters={'sampling_depth': Int % Range(1, None),
                'metric': Str % Choices(alpha_metrics['NONPHYLO']['IMPL'] |
                                        alpha_metrics['NONPHYLO']['UNIMPL'] |
                                        alpha_metrics['PHYLO']['IMPL'] |
                                        alpha_metrics['PHYLO']['UNIMPL']),
                'n': Int % Range(1, None),
                'average_method': Str % Choices(['median' , 'mean' , 'mode'])},
    outputs={},
    input_descriptions={},
    parameter_descriptions={},
    output_descriptions={},
    name={},
    description={}
)
