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

    def test_bootstrap_iteration_filters_samples(self):
        t = Table(np.array([[0, 1, 3], [1, 1, 2]]),
                  ['O1', 'O2'],
                  ['S1', 'S2', 'S3'])

        observed = _bootstrap_iteration(t, 6)
        self.assertTrue(observed.is_empty())

        observed = _bootstrap_iteration(t, 5)
        self.assertEqual(list(observed.ids(axis='sample')), ['S3'])

        observed = _bootstrap_iteration(t, 2)
        self.assertEqual(list(observed.ids(axis='sample')), ['S2', 'S3'])

        observed = _bootstrap_iteration(t, 1)
        self.assertEqual(list(observed.ids(axis='sample')), ['S1', 'S2', 'S3'])

    def test_bootstrap_iteration_obtains_expected_counts(self):
        t = Table(np.array([[0, 10, 30], [1, 10, 20]]),
                  ['O1', 'O2'],
                  ['S1', 'S2', 'S3'])

        observed = _bootstrap_iteration(t, 1)
        self.assertEqual(list(observed.sum(axis="sample")), [1., 1., 1.])

        observed = _bootstrap_iteration(t, 10)
        self.assertEqual(list(observed.sum(axis="sample")), [10., 10.])

        observed = _bootstrap_iteration(t, 19)
        self.assertEqual(list(observed.sum(axis="sample")), [19., 19.])

        observed = _bootstrap_iteration(t, 25)
        self.assertEqual(list(observed.sum(axis="sample")), [25.])

        observed = _bootstrap_iteration(t, 49)
        self.assertEqual(list(observed.sum(axis="sample")), [49.])


if __name__ == "__main__":
    main()
