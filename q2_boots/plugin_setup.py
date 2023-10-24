# ----------------------------------------------------------------------------
# Copyright (c) 2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from qiime2.plugin import (Plugin, Int, Range, Collection,
                           Citations)
from q2_types.feature_table import (
    FeatureTable, Frequency
)

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
    outputs={'subsampled_tables': FeatureTable[Frequency]},
    input_descriptions={},
    parameter_descriptions={},
    output_descriptions={},
    name='Bootstrap Iteration',
    description=''
)

plugin.pipelines.register_function(
    function=q2_boots.bootstrap,
    inputs={'table': FeatureTable[Frequency]},
    parameters={'sampling_depth': Int % Range(1, None),
                'n': Int % Range(1, None)},
    outputs={'bootstrapped_tables': Collection[FeatureTable[Frequency]]},
    input_descriptions={},
    parameter_descriptions={},
    output_descriptions={},
    name='Bootstrap',
    description=''
)
