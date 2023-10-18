# ----------------------------------------------------------------------------
# Copyright (c) 2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from qiime2.plugin import (Plugin, Int, Float, Range, Metadata, Str, Bool,
                           Choices, MetadataColumn, Categorical, List,
                           Citations, TypeMatch, TypeMap, Collection)

import q2_boots
Citations = Citations.load('citations.bub', package='q2_boots')
plugin = Plugin(
    name='boots',
    version=q2_boots.__version__
)