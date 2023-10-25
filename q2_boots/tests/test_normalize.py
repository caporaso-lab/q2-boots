# ----------------------------------------------------------------------------
# Copyright (c) 2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from unittest import TestCase, main

import numpy as np
from biom.table import Table

from q2_boots import _bootstrap_iteration


class TestBootstrapIteration(TestCase):

    def test_bootstrap_iteration(self):
        t = Table(np.array([[0, 1, 3], [1, 1, 2]]),
                  ['O1', 'O2'],
                  ['S1', 'S2', 'S3'])
        a = _bootstrap_iteration(t, 2)
        self.assertEqual(a.shape, (2, 2))


if __name__ == "__main__":
    main()
