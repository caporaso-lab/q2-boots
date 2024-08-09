# ----------------------------------------------------------------------------
# Copyright (c) 2024, Caporaso Lab (https://cap-lab.bio).
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import pandas as pd
import pandas.testing as pdt
from skbio import TreeNode

import qiime2
from qiime2.plugin.testing import TestPluginBase

from q2_boots import alpha_average


class AlphaAverageTests(TestPluginBase):
    package = 'q2_boots'

    def test_median(self):
        vector1 = pd.Series([1., 200.,], index=['S1', 'S2'], name='x')
        vector2 = pd.Series([3., 300.,], index=['S1', 'S2'], name='x')
        vector3 = pd.Series([900., 3000.,], index=['S1', 'S2'], name='x')
        vector_collection = dict(enumerate([vector1, vector2, vector3]))

        observed = alpha_average(vector_collection, average_method='median')
        expected = pd.Series([3., 300.,], index=['S1', 'S2'], name='x')
        pdt.assert_series_equal(observed, expected)

    def test_mean(self):
        vector1 = pd.Series([1., 200.,], index=['S1', 'S2'], name='x')
        vector2 = pd.Series([3., 300.,], index=['S1', 'S2'], name='x')
        vector3 = pd.Series([900., 3000.,], index=['S1', 'S2'], name='x')
        vector_collection = dict(enumerate([vector1, vector2, vector3]))

        observed = alpha_average(vector_collection, average_method='mean')
        expected = pd.Series([904./3, 3500./3,],
                             index=['S1', 'S2'],
                             name='x')
        pdt.assert_series_equal(observed, expected)

    def test_invalid_average_method(self):
        vector1 = pd.Series([1., 200.,], index=['S1', 'S2'], name='x')
        vector2 = pd.Series([3., 300.,], index=['S1', 'S2'], name='x')
        vector3 = pd.Series([900., 3000.,], index=['S1', 'S2'], name='x')
        vector_collection = dict(enumerate([vector1, vector2, vector3]))

        with self.assertRaisesRegex(KeyError, "'w'.*'median' and 'mean'"):
            alpha_average(vector_collection, average_method='w')


class AlphaCollectionTests(TestPluginBase):
    package = 'q2_boots'

    def setUp(self):
        super().setUp()
        self.alpha_collection_pipeline = \
            self.plugin.pipelines['alpha_collection']

    def test_alpha_collection_w_replacement(self):
        table1 = pd.DataFrame(data=[[1, 1], [0, 4]],
                              columns=['F1', 'F2'],
                              index=['S1', 'S2'])
        table1 = qiime2.Artifact.import_data(
            "FeatureTable[Frequency]", table1, view_type=pd.DataFrame
        )

        # at a sampling depth of 1, with table1 as input, there is one possible
        # outcome. confirm that we observe it.
        observed, = self.alpha_collection_pipeline(
            table=table1, sampling_depth=1, metric='observed_features', n=42,
            replacement=True)
        self.assertEqual(len(observed), 42)

        expected_series = pd.Series([1, 1],
                                    index=['S1', 'S2'],
                                    name='observed_features')
        for alpha_vector in observed.values():
            observed_series = alpha_vector.view(pd.Series)
            pdt.assert_series_equal(observed_series, expected_series)

        # at a sampling depth of 2, with table1 as input and sampling with
        # replacement, there are two possible outcomes. confirm that in 100
        # iterations each is observed at least once and no other outputs are
        # observed.
        observed, = self.alpha_collection_pipeline(
            table=table1, sampling_depth=2, metric='observed_features', n=100,
            replacement=True)
        self.assertEqual(len(observed), 100)

        possible_series1 = pd.Series([1, 1],
                                     index=['S1', 'S2'],
                                     name='observed_features')
        possible_series2 = pd.Series([2, 1],
                                     index=['S1', 'S2'],
                                     name='observed_features')
        count_possible_series1_observed = 0
        count_possible_series2_observed = 0
        other_series_observed = False
        for alpha_vector in observed.values():
            observed_series = alpha_vector.view(pd.Series)
            if observed_series.equals(possible_series1):
                count_possible_series1_observed += 1
            elif observed_series.equals(possible_series2):
                count_possible_series2_observed += 1
            else:
                other_series_observed = True
        self.assertTrue(count_possible_series1_observed > 0)
        self.assertTrue(count_possible_series2_observed > 0)
        self.assertFalse(other_series_observed)

    def test_alpha_collection_wo_replacement(self):
        table1 = pd.DataFrame(data=[[1, 1], [0, 4]],
                              columns=['F1', 'F2'],
                              index=['S1', 'S2'])
        table1 = qiime2.Artifact.import_data(
            "FeatureTable[Frequency]", table1, view_type=pd.DataFrame
        )

        # at a sampling depth of 2, with table1 as input and sampling without
        # replacement, there is one possible outcome. confirm that we always
        # observe it.
        observed, = self.alpha_collection_pipeline(
            table=table1, sampling_depth=2, metric='observed_features', n=88,
            replacement=False)
        self.assertEqual(len(observed), 88)

        expected_series = pd.Series([2, 1],
                                    index=['S1', 'S2'],
                                    name='observed_features')
        for alpha_vector in observed.values():
            observed_series = alpha_vector.view(pd.Series)
            pdt.assert_series_equal(observed_series, expected_series)

    def test_alpha_collection_phylogenetic(self):
        table1 = pd.DataFrame(data=[[1, 1], [0, 4]],
                              columns=['F1', 'F2'],
                              index=['S1', 'S2'])
        table1 = qiime2.Artifact.import_data(
            "FeatureTable[Frequency]", table1, view_type=pd.DataFrame
        )
        phylogeny = TreeNode.read(["((F1:1.0,F2:1.0):2.0);"])
        phylogeny = qiime2.Artifact.import_data(
            "Phylogeny[Rooted]", phylogeny)

        # at a sampling depth of 2, with table1 as input and sampling without
        # replacement, there is one possible outcome. confirm that we always
        # observe it.
        observed, = self.alpha_collection_pipeline(
            table=table1, sampling_depth=2, metric='faith_pd', n=88,
            replacement=False, phylogeny=phylogeny)
        self.assertEqual(len(observed), 88)

        expected_series = pd.Series([4, 3],
                                    index=['S1', 'S2'],
                                    name='faith_pd')
        for alpha_vector in observed.values():
            observed_series = alpha_vector.view(pd.Series)
            pdt.assert_series_equal(observed_series, expected_series)

    def test_alpha_collection_invalid_input(self):
        table1 = pd.DataFrame(data=[[1, 1], [0, 4]],
                              columns=['F1', 'F2'],
                              index=['S1', 'S2'])
        table1 = qiime2.Artifact.import_data(
            "FeatureTable[Frequency]", table1, view_type=pd.DataFrame
        )

        with self.assertRaisesRegex(ValueError, "_pd requires a phy"):
            self.alpha_collection_pipeline(table=table1, sampling_depth=2,
                                           metric='faith_pd', n=88,
                                           replacement=False)

        with self.assertRaisesRegex(TypeError, "xyz"):
            self.alpha_collection_pipeline(table=table1, sampling_depth=2,
                                           metric='xyz', n=88,
                                           replacement=False)


class AlphaTests(TestPluginBase):
    package = 'q2_boots'

    def setUp(self):
        super().setUp()
        self.alpha_pipeline = self.plugin.pipelines['alpha']

    def test_alpha_w_replacement(self):
        table1 = pd.DataFrame(data=[[1, 1], [0, 4]],
                              columns=['F1', 'F2'],
                              index=['S1', 'S2'])
        table1 = qiime2.Artifact.import_data(
            "FeatureTable[Frequency]", table1, view_type=pd.DataFrame
        )

        # at a sampling depth of 1, with table1 as input, there is one possible
        # outcome. confirm that we observe it.
        observed, = self.alpha_pipeline(
            table=table1, sampling_depth=1, metric='observed_features', n=42,
            replacement=True)
        observed_series = observed.view(pd.Series)

        expected_series = pd.Series([1., 1.],
                                    index=['S1', 'S2'],
                                    name='observed_features')
        pdt.assert_series_equal(observed_series, expected_series)

        # at a sampling depth of 2, with table1 as input, sampling with
        # replacement, and averaging with median, there are two possible
        # outcomes for S1 and one possible outcome for S2. confirm that in 99
        # iterations we observe one of the possible values for S1 and the
        # expected value for S2
        observed, = self.alpha_pipeline(
            table=table1, sampling_depth=2, metric='observed_features',
            n=99, replacement=True)
        observed_series = observed.view(pd.Series)
        # note that because n is an odd number, the median for S1 will always
        # be one of the actual values that were observed (as opposed to possibly
        # being 1.5, if n was an even number).
        self.assertTrue(observed_series['S1'] == 1.0 or
                        observed_series['S1'] == 2.0)
        self.assertEqual(observed_series['S2'], 1.0)

        # at a sampling depth of 2, with table1 as input, sampling with
        # replacement, and averaging with mean, there are many possible
        # outcomes for S1 and one possible outcome for S2. confirm that in
        # 100 iterations S1 is always in the expected range and S2 always has
        # the expected value
        observed, = self.alpha_pipeline(
            table=table1, sampling_depth=2, metric='observed_features',
            n=100, replacement=True, average_method='mean')
        observed_series = observed.view(pd.Series)
        # note that b/c these are <, not <=, we know that identical tables
        # were not always obtained from the resampling step
        self.assertTrue(1.0 < observed_series['S1'] < 2.0)
        self.assertEqual(observed_series['S2'], 1.0)

    def test_alpha_wo_replacement(self):
        table1 = pd.DataFrame(data=[[1, 1], [0, 4]],
                              columns=['F1', 'F2'],
                              index=['S1', 'S2'])
        table1 = qiime2.Artifact.import_data(
            "FeatureTable[Frequency]", table1, view_type=pd.DataFrame
        )

        # at a sampling depth of 1, with table1 as input, there is one possible
        # outcome. confirm that we observe it.
        observed, = self.alpha_pipeline(
            table=table1, sampling_depth=1, metric='observed_features', n=42,
            replacement=False)
        observed_series = observed.view(pd.Series)

        expected_series = pd.Series([1., 1.],
                                    index=['S1', 'S2'],
                                    name='observed_features')
        pdt.assert_series_equal(observed_series, expected_series)

        # at a sampling depth of 2, with table1 as input, sampling without
        # replacement, and averaging with median, there is one possible outcome.
        # confirm that we observe it.
        observed, = self.alpha_pipeline(
            table=table1, sampling_depth=2, metric='observed_features',
            n=100, replacement=False)
        observed_series = observed.view(pd.Series)

        expected_series = pd.Series([2., 1.],
                                    index=['S1', 'S2'],
                                    name='observed_features')
        pdt.assert_series_equal(observed_series, expected_series)

        # at a sampling depth of 2, with table1 as input, sampling without
        # replacement, and averaging with mean, there is one possible outcome.
        # confirm that we observe it.
        observed, = self.alpha_pipeline(
            table=table1, sampling_depth=2, metric='observed_features',
            n=100, replacement=False, average_method='mean')
        observed_series = observed.view(pd.Series)

        expected_series = pd.Series([2., 1.],
                                    index=['S1', 'S2'],
                                    name='observed_features')
        pdt.assert_series_equal(observed_series, expected_series)
