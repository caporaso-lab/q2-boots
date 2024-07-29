# ----------------------------------------------------------------------------
# Copyright (c) 2024, Caporaso Lab (https://cap-lab.bio).
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from io import StringIO

import numpy as np
import pandas as pd
import pandas.testing as pdt
import skbio
from biom.table import Table

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

        observed = alpha_average(vector_collection, average_method='median')
        expected = pd.Series([904./3, 3500./3,],
                             index=['S1', 'S2'],
                             name='x')
        pdt.assert_series_equal(observed, expected)

    def test_invalid_average_method(self):
        vector1 = pd.Series([1., 200.,], index=['S1', 'S2'], name='x')
        vector2 = pd.Series([3., 300.,], index=['S1', 'S2'], name='x')
        vector3 = pd.Series([900., 3000.,], index=['S1', 'S2'], name='x')
        vector_collection = dict(enumerate([vector1, vector2, vector3]))

        with self.assertRaisesRegex(KeyError, "median or mean"):
            alpha_average(vector_collection, average_method='w')


class AlphaCollectionTests(TestPluginBase):
    package = 'q2_boots'

    def setUp(self):
        super().setUp()
        self.alpha_pipeline = self.plugin.pipelines['alpha']
        self.alpha_collection_pipeline = \
            self.plugin.pipelines['alpha_collection']

    def test_alpha_w_replacement(self):
        table1 = pd.DataFrame(data=[[0, 2], [3, 2]],
                              columns=['F1', 'F2'],
                              index=['S1', 'S2'])
        table1 = qiime2.Artifact.import_data(
            "FeatureTable[Frequency]", table1, view_type=pd.DataFrame
        )

        observed = self.alpha_pipeline(
            table=table1, sampling_depth=1, metric='observed_features', n=10,
            replacement=True)
        self.assertEqual(len(observed), 1)
        observed_series = qiime2.Artifact.view(observed[0], view_type=pd.Series)
        expected_series = pd.Series([1., 1.],
                                    index=['S1', 'S2'],
                                    name='observed_features')
        pdt.assert_series_equal(observed_series, expected_series)

        # at a sampling depth of 2, with table1 as input, there are only two
        # possible outputs. each will occur with a probability of 0.5 on each
        # iteration. confirm that in 100 iterations, both are observed at least
        # once, and no other outputs are observed
        possible_series1 = pd.Series([1., 1.],
                                     index=['S1', 'S2'],
                                     name='observed_features')
        possible_series2 = pd.Series([1., 2.],
                                     index=['S1', 'S2'],
                                     name='observed_features')
        count_possible_series1_observed = 0
        count_possible_series2_observed = 0
        other_series_observed = False
        for _ in range(20):
            observed = self.alpha_pipeline(
                table=table1, sampling_depth=2, metric='observed_features',
                n=2, replacement=True)

            observed_series = qiime2.Artifact.view(
                observed[0], view_type=pd.Series)

            if observed_series.equals(possible_series1):
                count_possible_series1_observed += 1
            elif observed_series.equals(possible_series2):
                count_possible_series2_observed += 1
            else:
                other_series_observed = True

        self.assertTrue(count_possible_series1_observed > 0)
        self.assertTrue(count_possible_series2_observed > 0)
        self.assertFalse(other_series_observed)

        ## Hmm, the above passes, but I'm actually not sure why. If the default
        ## averaging method is median, I would expect this S2 to sometimes have
        ## a value of 1.5. Pick up here.
        self.assertTrue(False)

    def test_alpha_wo_replacement(self):
        table1 = pd.DataFrame(data=[[0, 1], [1, 1], [3, 2]],
                              columns=['F1', 'F2'],
                              index=['S1', 'S2', 'S3'])
        table1 = qiime2.Artifact.import_data(
            "FeatureTable[Frequency]", table1, view_type=pd.DataFrame
        )

        observed = self.alpha_pipeline(
            table=table1, sampling_depth=1, metric='observed_features', n=10,
            replacement=True)

        self.assertEqual(len(observed), 1)
        observed_series = qiime2.Artifact.view(observed[0], view_type=pd.Series)
        expected_series = pd.Series([1., 1., 1.],
                                    index=['S1', 'S2', 'S3'],
                                    name='observed_features')

        pdt.assert_series_equal(observed_series, expected_series)

    def test_range_non_phylo(self):
        t = Table(np.array([[0, 1, 3], [1, 1, 2], [1, 3, 2]]),
                  ['O1', 'O2', 'O3'],
                  ['S1', 'S2', 'S3'])
        t = qiime2.Artifact.import_data('FeatureTable[Frequency]', t)
        output, = self.alpha_pipeline(table=t, sampling_depth=1,
                             metric='shannon',
                             n=10, replacement=True)
        output: pd.Series = qiime2.Artifact.view(output, pd.Series)

        collection, = self.alpha_collection_pipeline(
            table=t, sampling_depth=1,
            metric='shannon', n=10, replacement=True
        )

        self.assertTrue(self._range_check(output, collection.values()))

    def test_range_phylo(self):
        with StringIO('(O1:0.3, O2:0.2, O3:0.1, O4:0.2)root;') as f:
            phylogeny = skbio.read(f, format='newick', into=skbio.TreeNode)

        phylogeny = qiime2.Artifact.import_data(
            'Phylogeny[Rooted]',
            phylogeny
        )
        t = Table(np.array([[0, 1, 3], [1, 1, 2], [1, 3, 2]]),
                  ['O1', 'O2', 'O3'],
                  ['S1', 'S2', 'S3'])
        t = qiime2.Artifact.import_data('FeatureTable[Frequency]', t)
        output, = self.alpha_pipeline(table=t, sampling_depth=1,
                             metric='pielou_e',
                             phylogeny=phylogeny,
                             n=10, replacement=True)
        output: pd.Series = qiime2.Artifact.view(output, pd.Series)

        collection, = self.alpha_collection_pipeline(
            table=t, sampling_depth=1,
            phylogeny=phylogeny,
            metric='pielou_e', n=10, replacement=True
        )

        self.assertTrue(self._range_check(output, collection.values()))

    def test_phylo_metric_no_phylo(self):
        t = Table(np.array([[0, 1, 3], [1, 1, 2], [1, 3, 2]]),
                  ['01', '02', '03'],
                  ['S1', 'S2', 'S3'])
        t = qiime2.Artifact.import_data('FeatureTable[Frequency]', t)

        with self.assertRaisesRegex(ValueError, 'non-phylogenic metric'):
            self.alpha_pipeline(table=t, sampling_depth=1,
                       metric='faith_pd',
                       n=10, replacement=True)

    def test_non_phylo_metric_with_phylo(self):
        with StringIO('(O1:0.3, O2:0.2, O3:0.1, O4:0.2)root;') as f:
            phylogeny = skbio.read(f, format='newick', into=skbio.TreeNode)

        phylogeny = qiime2.Artifact.import_data(
            'Phylogeny[Rooted]',
            phylogeny
        )
        t = Table(np.array([[0, 1, 3], [1, 1, 2], [1, 3, 2]]),
                  ['01', '02', '03'],
                  ['S1', 'S2', 'S3'])
        t = qiime2.Artifact.import_data('FeatureTable[Frequency]', t)
        # assert no error is thrown, and that phylogeny is just thrown out
        self.alpha_pipeline(table=t, sampling_depth=1,
                   metric='shannon',
                   n=10,
                   phylogeny=phylogeny, replacement=True)
        self.assertTrue(True)

    def _range_check(self, output, div_collection):

        df = None
        i = 0

        all_series = [qiime2.Artifact.view(x, pd.Series) for x in div_collection]

        for series in all_series:
            if df is None:
                df = pd.DataFrame(series)
            else:
                series.name = i
                df.join(series)
            i += 1

        max_value = df.max(axis=1)
        min_value = df.min(axis=1)

        for i in range(len(output)):
            if output[i] > max_value[i] or output[i] < min_value[i]:
                return False
        return True


class TestAlphaBootstrapRepresentative(TestPluginBase):

    package = 'q2_boots'

    def setUp(self):
        super().setUp()
        self.alpha_collection = self.plugin.pipelines[
            'alpha_collection']

    def check_samples(self, indices, data: pd.Series):
        for index in indices:
            self.assertTrue(index in data.index)

    def test_basic(self):
        t = Table(np.array([[0, 1, 3], [1, 1, 2], [1, 3, 2]]),
                  ['01', '02', '03'],
                  ['S1', 'S2', 'S3'])
        t = qiime2.Artifact.import_data('FeatureTable[Frequency]', t)
        output = self.alpha_collection(table=t, sampling_depth=1,
                                       metric='shannon',
                                       n=10, replacement=True)

        self.assertEqual(len(output[0]), 10)
        index = ['S1', 'S2', 'S3']

        for series in output[0].values():
            series = qiime2.Artifact.view(series, pd.Series)
            self.check_samples(index, series)
