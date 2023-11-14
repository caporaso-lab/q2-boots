# ----------------------------------------------------------------------------
# Copyright (c) 2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from qiime2.plugin import (Plugin, Int, Range, Collection,
                           Citations, Str, Choices)
from q2_types.feature_table import (
    FeatureTable, Frequency
)
from q2_types.sample_data import AlphaDiversity, SampleData

from q2_types.tree import Phylogeny, Rooted

from q2_diversity_lib.alpha import METRICS

import q2_boots
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
    inputs={'table': FeatureTable[Frequency]},
    parameters={'sampling_depth': Int % Range(1, None)},
    outputs={'subsampled_table': FeatureTable[Frequency]},
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
    outputs={'subsampled_tables': Collection[FeatureTable[Frequency]]},
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
    function=q2_boots.alpha_bootstrap,
    inputs={'table': FeatureTable[Frequency],
            'phylogeny': Phylogeny[Rooted]},
    parameters={'sampling_depth': Int % Range(1, None),
                'metric': Str % Choices(METRICS['NONPHYLO']['IMPL'] |
                                        METRICS['NONPHYLO']['UNIMPL']),
                'n': Int % Range(1, None)},
    outputs={'sample_data': Collection[SampleData[AlphaDiversity]]},
    input_descriptions={'table': 'The table to be diversified',
                        'phylogeny': ''},
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
    description=''
)

plugin.pipelines.register_function(
    function=q2_boots.alpha_bootstrap_representative,
    inputs={'table': FeatureTable[Frequency],
            'phylogeny': Phylogeny[Rooted]},
    parameters={'sampling_depth': Int % Range(1, None),
                'metric': Str % Choices(METRICS['NONPHYLO']['IMPL'] |
                                        METRICS['NONPHYLO']['UNIMPL']),
                'n': Int % Range(1, None),
                'average_method': Str % Choices('median' | 'mean' | 'mode')},
    outputs={'sample_data': SampleData[AlphaDiversity]},
    input_descriptions={'table': 'The table to be diversified',
                        'phylogeny': ''},
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
