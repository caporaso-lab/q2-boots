# ----------------------------------------------------------------------------
# Copyright (c) 2024, Caporaso Lab (https://cap-lab.bio).
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
from q2_boots._beta import per_cell_average, get_medoid
import numpy as np
from skbio import DistanceMatrix
import skbio
from io import StringIO


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
        self.beta_collection = self.plugin.pipelines['beta_collection']

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

    def test_range_non_phylo(self):
        t = Table(np.array([[3000, 4000, 4151],
                            [1611, 3016, 2313], [3761, 2861, 2091]]),
                  ['01', '02', '03'],
                  ['S1', 'S2', 'S3'])
        t = Artifact.import_data('FeatureTable[Frequency]', t)
        output = self.beta(table=t,
                           metric='jaccard',
                           sampling_depth=500,
                           n=10,
                           representative='medoid'
                           )

        self.assertEqual(len(output), 1)

        output: DistanceMatrix = Artifact.view(output[0], DistanceMatrix)

        self.assertTrue('S1' in output.ids)
        self.assertTrue('S2' in output.ids)
        self.assertTrue('S3' in output.ids)

        output = output.to_data_frame()

        collection = self.beta_collection(
            table=t,
            metric='jaccard',
            sampling_depth=500,
            n=10,
        )

        index = output.index
        columns = output.columns

        collection = collection[0].values()
        mins, maxes = self.get_mins_and_maxes(collection)

        for col in columns:
            for row in index:
                self.assertTrue(
                    output[row][col] >= mins[row][col] and
                    output[row][col] <= maxes[row][col]
                )

    def test_range_phylo(self):
        t = Table(np.array([[3000, 4000, 4151],
                            [1611, 3016, 2313], [3761, 2861, 2091]]),
                  ['O1', 'O2', 'O3'],
                  ['S1', 'S2', 'S3'])
        t = Artifact.import_data('FeatureTable[Frequency]', t)
        with StringIO('(O1:0.3, O2:0.2, O3:0.1, O4:0.2)root;') as f:
            phylogeny = skbio.read(f, format='newick', into=skbio.TreeNode)

        phylogeny = Artifact.import_data(
            'Phylogeny[Rooted]',
            phylogeny
        )
        output = self.beta(table=t,
                           metric='weighted_unifrac',
                           sampling_depth=500,
                           n=10,
                           representative='medoid',
                           phylogeny=phylogeny
                           )

        self.assertEqual(len(output), 1)

        output: DistanceMatrix = Artifact.view(output[0], DistanceMatrix)

        self.assertTrue('S1' in output.ids)
        self.assertTrue('S2' in output.ids)
        self.assertTrue('S3' in output.ids)

        output = output.to_data_frame()

        collection = self.beta_collection(
            table=t,
            metric='weighted_unifrac',
            sampling_depth=500,
            n=10,
            phylogeny=phylogeny,
        )

        index = output.index
        columns = output.columns

        collection = collection[0].values()
        mins, maxes = self.get_mins_and_maxes(collection)

        for col in columns:
            for row in index:
                self.assertTrue(
                    output[row][col] >= mins[row][col] and
                    output[row][col] <= maxes[row][col]
                )

    def test_phylo_metric_no_phylo(self):
        t = Table(np.array([[0, 1, 3], [1, 0, 2], [1, 3, 0]]),
                  ['01', '02', '03'],
                  ['S1', 'S2', 'S3'])
        t = Artifact.import_data('FeatureTable[Frequency]', t)

        with self.assertRaisesRegex(ValueError, 'non-phylogenic metric'):
            self.beta(table=t,
                      metric='weighted_unifrac',
                      sampling_depth=1,
                      n=10,
                      representative='medoid'
                      )

    def test_non_phylo_metric_phylo(self):
        t = Table(np.array([[3000, 4000, 4151],
                            [1611, 3016, 2313], [3761, 2861, 2091]]),
                  ['O1', 'O2', 'O3'],
                  ['S1', 'S2', 'S3'])
        t = Artifact.import_data('FeatureTable[Frequency]', t)
        with StringIO('(O1:0.3, O2:0.2, O3:0.1, O4:0.2)root;') as f:
            phylogeny = skbio.read(f, format='newick', into=skbio.TreeNode)

        phylogeny = Artifact.import_data(
            'Phylogeny[Rooted]',
            phylogeny
        )
        # just assert no value is raised
        self.beta(table=t,
                  metric='jaccard',
                  sampling_depth=5,
                  n=10,
                  representative='medoid',
                  phylogeny=phylogeny
                  )
        self.assertTrue(True)

    def get_mins_and_maxes(self, dms):

        dfs = []
        for dm in dms:
            dm = Artifact.view(dm, DistanceMatrix)
            dfs.append(dm.to_data_frame())
        x, y = dfs[0].shape

        columns = dfs[0].columns
        index = dfs[0].index

        mins = pd.DataFrame(index=index, columns=columns)
        maxes = pd.DataFrame(index=index, columns=columns)

        for col in columns:
            for row in index:
                vals = []
                for df in dfs:
                    vals.append(df[row][col])
                mins[row][col] = min(vals)
                maxes[row][col] = max(vals)

        return mins, maxes


if __name__ == "__main__":
    main()
