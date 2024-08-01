# ----------------------------------------------------------------------------
# Copyright (c) 2024, Caporaso Lab (https://cap-lab.bio).
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import qiime2
from qiime2.plugin.testing import TestPluginBase
import pandas as pd
import pandas.testing as pdt
import skbio


class CoreMetricsTests(TestPluginBase):

    package = 'q2_boots'

    def setUp(self):
        super().setUp()
        self.core_metrics = self.plugin.pipelines['core_metrics']
        table1 = pd.DataFrame(data=[[1, 1], [0, 4]],
                              columns=['F1', 'F2'],
                              index=['S1', 'S2'])
        self.table1 = qiime2.Artifact.import_data(
            "FeatureTable[Frequency]", table1, view_type=pd.DataFrame)

        phylogeny = skbio.TreeNode.read(["((F1:1.0,F2:1.0):2.0);"])
        self.phylogeny = qiime2.Artifact.import_data(
            "Phylogeny[Rooted]", phylogeny)

        metadata = pd.DataFrame(['not', 'of', 'interest'],
                                index=['S1', 'S2', 'S3'],
                                columns=['blank'])
        metadata.index.name = 'sample-id'
        self.metadata = qiime2.Metadata(metadata)

    def test_core_metrics_wo_replacement(self):
        output = self.core_metrics(table=self.table1,
                                   sampling_depth=2,
                                   metadata=self.metadata,
                                   n_jobs=1,
                                   replacement=False,
                                   n=10)
        # n tables returned
        self.assertEqual(len(output[0]), 10)
        expected_table = pd.DataFrame(data=[[1.0, 1.0], [0.0, 2.0]],
                                      columns=['F1', 'F2'],
                                      index=['S1', 'S2'])
        for e in output[0].values():
            observed_table = qiime2.Artifact.view(e, view_type=pd.DataFrame)
            pdt.assert_frame_equal(observed_table, expected_table)

        # expected alpha vectors returned
        self.assertEqual(set(output[1].keys()),
                         set(['observed_features', 'pielou_evenness',
                              'shannon_entropy']))
        expected_obs_features = pd.Series([2.0, 1.0],
                                          index=['S1', 'S2'],
                                          name='observed_features')
        observed_obs_features = qiime2.Artifact.view(
            output[1]['observed_features'], view_type=pd.Series)
        pdt.assert_series_equal(observed_obs_features, expected_obs_features)

        # expected dms, pcoas, and plots returned
        self.assertEqual(set(output[2].keys()), set(['jaccard', 'braycurtis']))
        self.assertEqual(set(output[3].keys()), set(['jaccard', 'braycurtis']))
        self.assertEqual(set(output[4].keys()), set(['jaccard', 'braycurtis']))
        expected_jaccard = skbio.DistanceMatrix([[0, 0.5], [0.5, 0]],
                                                ids=['S1', 'S2'])
        observed_jaccard = qiime2.Artifact.view(output[2]['jaccard'],
                                                view_type=skbio.DistanceMatrix)
        self.assertEqual(observed_jaccard, expected_jaccard)

    def test_core_metrics_w_replacement(self):
        output = self.core_metrics(table=self.table1,
                                   sampling_depth=2,
                                   metadata=self.metadata,
                                   n_jobs=1,
                                   replacement=True,
                                   n=100)
        # n tables returned
        self.assertEqual(len(output[0]), 100)
        possible_table1 = pd.DataFrame(data=[[1.0, 1.0], [0.0, 2.0]],
                                       columns=['F1', 'F2'],
                                       index=['S1', 'S2'])
        possible_table2 = pd.DataFrame(data=[[2.0, 0.0], [0.0, 2.0]],
                                       columns=['F1', 'F2'],
                                       index=['S1', 'S2'])
        possible_table3 = pd.DataFrame(data=[[2.0], [2.0]],
                                       columns=['F2'],
                                       index=['S1', 'S2'])
        count_possible_table1_observed = 0
        count_possible_table2_observed = 0
        count_possible_table3_observed = 0
        count_other_table_observed = 0
        for e in output[0].values():
            observed_table = qiime2.Artifact.view(e, view_type=pd.DataFrame)
            if observed_table.equals(possible_table1):
                count_possible_table1_observed += 1
            elif observed_table.equals(possible_table2):
                count_possible_table2_observed += 1
            elif observed_table.equals(possible_table3):
                count_possible_table3_observed += 1
            else:
                count_other_table_observed += 1
        self.assertTrue(count_possible_table1_observed > 0)
        self.assertTrue(count_possible_table2_observed > 0)
        self.assertTrue(count_possible_table3_observed > 0)
        self.assertEqual(count_other_table_observed, 0)

        # expected alpha vectors returned
        self.assertEqual(set(output[1].keys()),
                         set(['observed_features', 'pielou_evenness',
                              'shannon_entropy']))
        observed_obs_features = qiime2.Artifact.view(
            output[1]['observed_features'], view_type=pd.Series)
        self.assertTrue(observed_obs_features['S1'] == 1.0 or
                        observed_obs_features['S1'] == 2.0)
        self.assertEqual(observed_obs_features['S2'], 1.0)

        # expected dms, pcoas, and plots returned
        self.assertEqual(set(output[2].keys()), set(['jaccard', 'braycurtis']))
        self.assertEqual(set(output[3].keys()), set(['jaccard', 'braycurtis']))
        self.assertEqual(set(output[4].keys()), set(['jaccard', 'braycurtis']))
        observed_jaccard = qiime2.Artifact.view(output[2]['jaccard'],
                                                view_type=skbio.DistanceMatrix)
        self.assertTrue(0.0 < observed_jaccard[('S1', 'S2')] < 1.0)

    def test_core_metrics_phylogenetic(self):
        output = self.core_metrics(table=self.table1,
                                   phylogeny=self.phylogeny,
                                   sampling_depth=2,
                                   metadata=self.metadata,
                                   n_jobs=1,
                                   replacement=False,
                                   n=10)
        # n tables returned
        self.assertEqual(len(output[0]), 10)
        expected_table = pd.DataFrame(data=[[1.0, 1.0], [0.0, 2.0]],
                                      columns=['F1', 'F2'],
                                      index=['S1', 'S2'])
        for e in output[0].values():
            observed_table = qiime2.Artifact.view(e, view_type=pd.DataFrame)
            pdt.assert_frame_equal(observed_table, expected_table)

        # expected alpha vectors returned
        self.assertEqual(set(output[1].keys()),
                         set(['observed_features', 'pielou_evenness',
                              'shannon_entropy', 'faith_pd']))
        expected_obs_features = pd.Series([2.0, 1.0],
                                          index=['S1', 'S2'],
                                          name='observed_features')
        observed_obs_features = qiime2.Artifact.view(
            output[1]['observed_features'], view_type=pd.Series)
        pdt.assert_series_equal(observed_obs_features, expected_obs_features)

        expected_faith_pd = pd.Series([4.0, 3.0],
                                      index=['S1', 'S2'],
                                      name='faith_pd')
        observed_faith_pd = qiime2.Artifact.view(
            output[1]['faith_pd'], view_type=pd.Series)
        pdt.assert_series_equal(observed_faith_pd, expected_faith_pd)

        # expected dms, pcoas, and plots returned
        expected_bdiv_keys = set(['jaccard', 'braycurtis',
                                  'unweighted_unifrac', 'weighted_unifrac'])
        self.assertEqual(set(output[2].keys()), expected_bdiv_keys)
        self.assertEqual(set(output[3].keys()), expected_bdiv_keys)
        self.assertEqual(set(output[4].keys()), expected_bdiv_keys)
        expected_jaccard = skbio.DistanceMatrix([[0, 0.5], [0.5, 0]],
                                                ids=['S1', 'S2'])
        observed_jaccard = qiime2.Artifact.view(output[2]['jaccard'],
                                                view_type=skbio.DistanceMatrix)
        self.assertEqual(observed_jaccard, expected_jaccard)

        expected_unweighted_unifrac = skbio.DistanceMatrix(
            [[0, 0.25], [0.25, 0]], ids=['S1', 'S2'])
        observed_unweighted_unifrac = qiime2.Artifact.view(
            output[2]['unweighted_unifrac'], view_type=skbio.DistanceMatrix)
        self.assertEqual(observed_unweighted_unifrac,
                         expected_unweighted_unifrac)
