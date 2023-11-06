# ----------------------------------------------------------------------------
# Copyright (c) 2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from ._normalize import _bootstrap_iteration, bootstrap
from .alpha import alpha_bootstrap, alpha_phylogenetic_bootstrap
from ._version import get_versions

__all__ = ['_bootstrap_iteration', 'bootstrap',
           'alpha_bootstrap', 'alpha_phylogenetic_bootstrap']

__version__ = get_versions()['version']
del get_versions
