# ----------------------------------------------------------------------------
# Copyright (c) 2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from ._normalize import _bootstrap_iteration, resample
from .alpha import (alpha,
                    alpha_collection,
                    alpha_average)
from ._version import get_versions
from .beta import (beta, get_medoid, beta_collection, per_cell_average, beta_average)
from .core_metrics import (core_metrics)

__all__ = ['_bootstrap_iteration', 'resample',
           'alpha',
           'alpha_collection',
           'beta',
           'core_metrics_alpha_bootstrap',
           'core_metrics',
           'get_medoid',
           'per_cell_average',
           'beta_collection',
           'beta_average',
           'alpha_average']

__version__ = get_versions()['version']
del get_versions

from . import _version
__version__ = _version.get_versions()['version']
