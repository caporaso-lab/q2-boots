# ----------------------------------------------------------------------------
# Copyright (c) 2024, Caporaso Lab (https://cap-lab.bio).
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from unittest import TestCase, main

import pandas as pd
import qiime2

import skbio

from qiime2.plugin.testing import TestPluginBase

from q2_boots import beta_average
from q2_boots._beta import _per_cell_average, _medoid


class BetaAverageTests(TestCase):

    def setUp(self):
        super().setUp()
        self.a = skbio.DistanceMatrix([[0, 2, 99],
                                       [2, 0, 1],
                                       [99, 1, 0]], ids=('S1', 'S2', 'S3'))
        self.b = skbio.DistanceMatrix([[0, 4, 1],
                                       [4, 0, 2],
                                       [1, 2, 0]], ids=('S1', 'S2', 'S3'))
        self.c = skbio.DistanceMatrix([[0, 6, 2],
                                       [6, 0, 3],
                                       [2, 3, 0]], ids=('S1', 'S2', 'S3'))
        self.dms = {'a': self.a, 'b': self.b, 'c': self.c}

    def test_non_metric_median(self):
        observed = beta_average(self.dms, "non-metric-median")
        exp = skbio.DistanceMatrix([[0.0, 4.0, 2.0],
                                    [4.0, 0.0, 2.0],
                                    [2.0, 2.0, 0.0]], ids=('S1', 'S2', 'S3'))

        self.assertEqual(observed, exp)

        self.assertEqual(beta_average({'a': self.a}, "non-metric-median"),
                         self.a)
        self.assertEqual(beta_average({'a1': self.a,
                                       'a2': self.a,
                                       'a3': self.a,
                                       'a4': self.a}, "non-metric-median"),
                         self.a)

    def test_non_metric_mean(self):
        observed = beta_average(self.dms, "non-metric-mean")
        exp = skbio.DistanceMatrix([[0.0, 4.0, 34.0],
                                    [4.0, 0.0, 2.0],
                                    [34.0, 2.0, 0.0]], ids=('S1', 'S2', 'S3'))

        self.assertEqual(observed, exp)

        self.assertEqual(beta_average({'a': self.a}, "non-metric-mean"),
                         self.a)
        self.assertEqual(beta_average({'a1': self.a,
                                       'a2': self.a,
                                       'a3': self.a,
                                       'a4': self.a}, "non-metric-mean"),
                         self.a)

    def test_medoid(self):
        observed = beta_average(self.dms, "medoid")
        exp = skbio.DistanceMatrix([[0, 6, 2],
                                    [6, 0, 3],
                                    [2, 3, 0]], ids=('S1', 'S2', 'S3'))

        self.assertEqual(observed, exp)

        self.assertEqual(beta_average({'a': self.a}, "medoid"),
                         self.a)
        self.assertEqual(beta_average({'a1': self.a,
                                       'a2': self.a,
                                       'a3': self.a,
                                       'a4': self.a}, "medoid"),
                         self.a)

    def test_invalid(self):
        with self.assertRaisesRegex(ValueError, "Unknown average method.*xyz"):
            beta_average(self.dms, "xyz")


class BetaAverageHelperTests(TestCase):

    def setUp(self):
        super().setUp()
        self.a = skbio.DistanceMatrix([[0, 2, 99],
                                       [2, 0, 1],
                                       [99, 1, 0]], ids=('S1', 'S2', 'S3'))
        self.b = skbio.DistanceMatrix([[0, 4, 1],
                                       [4, 0, 2],
                                       [1, 2, 0]], ids=('S1', 'S2', 'S3'))
        self.c = skbio.DistanceMatrix([[0, 6, 2],
                                       [6, 0, 3],
                                       [2, 3, 0]], ids=('S1', 'S2', 'S3'))
        self.dms = [self.a, self.b, self.c]

    def test_per_cell_median(self):
        observed = _per_cell_average(self.dms, 'median')
        exp = skbio.DistanceMatrix([[0.0, 4.0, 2.0],
                                    [4.0, 0.0, 2.0],
                                    [2.0, 2.0, 0.0]], ids=('S1', 'S2', 'S3'))

        self.assertEqual(observed, exp)

    def test_per_cell_mean(self):
        observed = _per_cell_average(self.dms, 'mean')
        exp = skbio.DistanceMatrix([[0.0, 4.0, 34.0],
                                    [4.0, 0.0, 2.0],
                                    [34.0, 2.0, 0.0]], ids=('S1', 'S2', 'S3'))

        self.assertEqual(observed, exp)

    def test_medoid(self):
        observed = _medoid(self.dms)
        self.assertEqual(observed, self.c)


class BetaCollectionTests(TestPluginBase):

    package = 'q2_boots'

    def setUp(self):
        super().setUp()
        self.beta_collection_pipeline = self.plugin.pipelines['beta_collection']

        table1 = pd.DataFrame(data=[[1, 1], [0, 4]],
                              columns=['F1', 'F2'],
                              index=['S1', 'S2'])
        self.table1 = qiime2.Artifact.import_data(
            "FeatureTable[Frequency]", table1, view_type=pd.DataFrame
        )

    def test_beta_collection_invalid_input(self):
        with self.assertRaisesRegex(ValueError, 'requires a phylogenetic tree'):
            self.beta_collection_pipeline(
                table=self.table1, metric='weighted_unifrac', sampling_depth=1,
                n=10, replacement=False)

    def test_beta_collection_w_replacement(self):
        # At a sampling depth of 2, with self.table1, and when sampling with
        # replacement, there are three possible Jaccard distance matrices.
        # Confirm that we see each of these at least once and no others.
        observed, = self.beta_collection_pipeline(
            table=self.table1, metric='jaccard', sampling_depth=2, n=100,
            replacement=True)
        self.assertEqual(len(observed), 100)

        possible_dm1 = skbio.DistanceMatrix([[0, 0.5], [0.5, 0]],
                                            ids=['S1', 'S2'])
        possible_dm2 = skbio.DistanceMatrix([[0, 1.0], [1.0, 0]],
                                            ids=['S1', 'S2'])
        possible_dm3 = skbio.DistanceMatrix([[0, 0.0], [0.0, 0]],
                                            ids=['S1', 'S2'])
        count_possible_dm1_observed = 0
        count_possible_dm2_observed = 0
        count_possible_dm3_observed = 0
        count_other_dm_observed = 0

        for o in observed.values():
            o = o.view(skbio.DistanceMatrix)
            if o == possible_dm1:
                count_possible_dm1_observed += 1
            elif o == possible_dm2:
                count_possible_dm2_observed += 1
            elif o == possible_dm3:
                count_possible_dm3_observed += 1
            else:
                count_other_dm_observed += 1
        self.assertTrue(count_possible_dm1_observed > 0)
        self.assertTrue(count_possible_dm2_observed > 0)
        self.assertTrue(count_possible_dm3_observed > 0)
        self.assertEqual(count_other_dm_observed, 0)

    def test_beta_collection_wo_replacement(self):
        # At a sampling depth of 2, with self.table1, and when sampling without
        # replacement, there is only one possible Jaccard distance matrix.
        # Confirm that we only ever see that distance matrix.
        expected = skbio.DistanceMatrix([[0, 0.5], [0.5, 0]], ids=['S1', 'S2'])

        observed, = self.beta_collection_pipeline(
            table=self.table1, metric='jaccard', sampling_depth=2, n=10,
            replacement=False)
        self.assertEqual(len(observed), 10)
        for o in observed.values():
            o = o.view(skbio.DistanceMatrix)
            self.assertEqual(o, expected)

    def test_beta_collection_phylogenetic(self):
        # At a sampling depth of 2, with self.table1, and when sampling without
        # replacement, there is only one possible unweighted UniFrac distance
        # matrix. Confirm that we only ever see that distance matrix.
        phylogeny = skbio.TreeNode.read(["((F1:1.0,F2:1.0):2.0);"])
        phylogeny = qiime2.Artifact.import_data(
            "Phylogeny[Rooted]", phylogeny)
        expected = skbio.DistanceMatrix([[0, 0.25], [0.25, 0]],
                                        ids=['S1', 'S2'])

        observed, = self.beta_collection_pipeline(
            table=self.table1, metric='unweighted_unifrac', phylogeny=phylogeny,
            sampling_depth=2, n=10, replacement=False)
        self.assertEqual(len(observed), 10)
        for o in observed.values():
            o = o.view(skbio.DistanceMatrix)
            self.assertEqual(o, expected)


class BetaTests(TestPluginBase):

    package = 'q2_boots'

    def setUp(self):
        super().setUp()
        self.beta_pipeline = self.plugin.pipelines['beta']

        table1 = pd.DataFrame(data=[[1, 1], [0, 4]],
                              columns=['F1', 'F2'],
                              index=['S1', 'S2'])
        self.table1 = qiime2.Artifact.import_data(
            "FeatureTable[Frequency]", table1, view_type=pd.DataFrame)

    def test_beta_w_replacement(self):
        # At a sampling depth of 2, with self.table1, and when sampling with
        # replacement, there are three possible Jaccard distance matrices.
        # Confirm the average is in range with all averaging methods.
        observed, = self.beta_pipeline(table=self.table1,
                                       metric='jaccard',
                                       sampling_depth=2,
                                       n=10,
                                       average_method='medoid',
                                       replacement=True)
        observed = observed.view(skbio.DistanceMatrix)
        self.assertTrue(observed[('S1', 'S2')] in [0.0, 0.5, 1.0],
                        msg=(f"Medoid value ({observed[('S1', 'S2')]}) is "
                             "not equal to 0.0, 0.5, 1.0."))

        observed, = self.beta_pipeline(table=self.table1,
                                       metric='jaccard',
                                       sampling_depth=2,
                                       n=9,
                                       average_method='non-metric-median',
                                       replacement=True)
        observed = observed.view(skbio.DistanceMatrix)
        # because n is odd, we should always observe an actual distance
        # between S1 and S2 as the median (as opposed to the mean of two
        # actual values)
        self.assertTrue(observed[('S1', 'S2')] in [0.0, 0.5, 1.0],
                        msg=(f"Median value ({observed[('S1', 'S2')]}) is "
                             "not equal to 0.0, 0.5 or 1.0."))

        observed, = self.beta_pipeline(table=self.table1,
                                       metric='jaccard',
                                       sampling_depth=2,
                                       n=20,
                                       average_method='non-metric-mean',
                                       replacement=True)
        observed = observed.view(skbio.DistanceMatrix)
        # This can occasionally fail, but it should be very infrequent
        # (e.g., if all 0.0s or 1.0s were observed as the distances)
        self.assertTrue(0.0 < observed[('S1', 'S2')] < 1.0,
                        msg=(f"Mean value ({observed[('S1', 'S2')]}) is "
                             "not in range to (0.0, 1.0)."))

    def test_beta_wo_replacement(self):
        # At a sampling depth of 2, with self.table1, and when sampling without
        # replacement, there is only one possible Jaccard distance matrix.
        # Confirm that we see it with all averaging methods.
        expected = skbio.DistanceMatrix([[0, 0.5], [0.5, 0]], ids=['S1', 'S2'])

        observed, = self.beta_pipeline(table=self.table1,
                                       metric='jaccard',
                                       sampling_depth=2,
                                       n=10,
                                       average_method='medoid',
                                       replacement=False)
        observed = observed.view(skbio.DistanceMatrix)
        self.assertEqual(observed, expected)

        observed, = self.beta_pipeline(table=self.table1,
                                       metric='jaccard',
                                       sampling_depth=2,
                                       n=10,
                                       average_method='non-metric-median',
                                       replacement=False)
        observed = observed.view(skbio.DistanceMatrix)
        self.assertEqual(observed, expected)

        observed, = self.beta_pipeline(table=self.table1,
                                       metric='jaccard',
                                       sampling_depth=2,
                                       n=10,
                                       average_method='non-metric-mean',
                                       replacement=False)
        observed = observed.view(skbio.DistanceMatrix)
        self.assertEqual(observed, expected)

    def test_invalid(self):
        with self.assertRaisesRegex(ValueError, 'requires a phylogenetic tree'):
            self.beta_pipeline(table=self.table1,
                               metric='weighted_unifrac',
                               sampling_depth=1,
                               n=10,
                               average_method='medoid',
                               replacement=False)

        with self.assertRaisesRegex(TypeError, 'xyz'):
            self.beta_pipeline(table=self.table1,
                               metric='jaccard',
                               sampling_depth=1,
                               n=10,
                               average_method='xyz',
                               replacement=False)


if __name__ == "__main__":
    main()
