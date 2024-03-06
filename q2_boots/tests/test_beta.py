# ----------------------------------------------------------------------------
# Copyright (c) 2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from unittest import TestCase, main
from qiime2.plugin.testing import TestPluginBase
import pandas as pd
from biom.table import Table
from qiime2 import Artifact
from q2_boots.beta import per_cell_average, get_medoid
import numpy as np
from skbio import DistanceMatrix


class TestAveraging(TestCase):

    def test_per_cell_median(self):

        a = pd.DataFrame([[1, 1, 1],
                          [1, 1, 1],
                          [1, 1, 1]])
        b = pd.DataFrame([[2, 2, 2],
                          [2, 2, 2],
                          [2, 2, 2]])
        c = pd.DataFrame([[3, 3, 3],
                          [3, 3, 3],
                          [3, 3, 3]])

        result = per_cell_average([a, b, c], 'median')
        exp = pd.DataFrame([[2.0, 2.0, 2.0],
                            [2.0, 2.0, 2.0],
                            [2.0, 2.0, 2.0]])

        pd.testing.assert_frame_equal(exp, result)

    def test_per_cell_mean(self):
        a = pd.DataFrame([[1, 0, 1],
                          [1, 0, 1],
                          [1, 0, 1]])
        b = pd.DataFrame([[2, 2, 2],
                          [2, 2, 2],
                          [2, 2, 2]])
        c = pd.DataFrame([[3, 1, 3],
                          [3, 1, 3],
                          [3, 1, 3]])

        result = per_cell_average([a, b, c], 'mean')
        exp = pd.DataFrame([[2.0, 1.0, 2.0],
                            [2.0, 1.0, 2.0],
                            [2.0, 1.0, 2.0]])

        pd.testing.assert_frame_equal(exp, result)

    def test_medoid(self):
        a = pd.DataFrame([[0, 0, 1],
                          [1, 0, 1],
                          [1, 0, 0]])
        b = pd.DataFrame([[0, 2, 2],
                          [2, 0, 2],
                          [2, 2, 0]])
        c = pd.DataFrame([[0, 1, 3],
                          [3, 0, 3],
                          [3, 1, 0]])

        result = get_medoid([a, b, c])
        pd.testing.assert_frame_equal(result, b)


class TestBeta(TestPluginBase):

    package = 'q2_boots'

    def setUp(self):
        super().setUp()
        self.beta = self.plugin.pipelines['beta']

    def test_basic(self):
        t = Table(np.array([[0, 1, 3], [1, 0, 2], [1, 3, 0]]),
                  ['01', '02', '03'],
                  ['S1', 'S2', 'S3'])
        t = Artifact.import_data('FeatureTable[Frequency]', t)
        output = self.beta(table=t,
                           metric='jaccard',
                           sampling_depth=1,
                           n=10,
                           representative='medoid'
                           )

        self.assertEqual(len(output), 1)

        output: DistanceMatrix = Artifact.view(output[0], DistanceMatrix)

        self.assertTrue('S1' in output.ids)
        self.assertTrue('S2' in output.ids)
        self.assertTrue('S3' in output.ids)

        output = output.to_data_frame()

        self.assertEqual(output.shape, (3, 3))


if __name__ == "__main__":
    main()
