# flake8: noqa
# ----------------------------------------------------------------------------
# Copyright (c) 2024, Caporaso Lab (https://cap-lab.bio).
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from ._resample import resample
from ._alpha import alpha, alpha_collection, alpha_average
from ._version import get_versions
from ._beta import beta, beta_collection, beta_average
from ._core_metrics import core_metrics
from . import _version

__all__ = ['resample',
           'alpha_average',
           'alpha_collection',
           'alpha',
           'beta_average'
           'beta_collection',
           'beta',
           'core_metrics']

__version__ = get_versions()['version']
del get_versions

__version__ = _version.get_versions()['version']
