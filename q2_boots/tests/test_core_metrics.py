# ----------------------------------------------------------------------------
# Copyright (c) 2024, Caporaso Lab (https://cap-lab.bio).
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import qiime2
from qiime2.plugin.testing import TestPluginBase
from qiime2.plugin import Visualization
import pandas as pd
import numpy.testing as npt
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
                                   replacement=False,
                                   n=10)
        # n tables returned
        self.assertEqual(len(output[0]), 10)
        expected_table = pd.DataFrame(data=[[1.0, 1.0], [0.0, 2.0]],
                                      columns=['F1', 'F2'],
                                      index=['S1', 'S2'])
        for e in output[0].values():
            observed_table = e.view(pd.DataFrame)
            pdt.assert_frame_equal(observed_table, expected_table)

        # expected alpha vectors returned
        skbio_lt_060_alpha_keys = set(
            ['observed_features', 'pielou_evenness', 'shannon_entropy'])
        skbio_gte_060_alpha_keys = set(
            ['observed_features', 'pielou_e', 'shannon'])
        self.assertTrue(set(output[1].keys()) == skbio_lt_060_alpha_keys or
                        set(output[1].keys()) == skbio_gte_060_alpha_keys)
        expected_obs_features = pd.Series([2.0, 1.0],
                                          index=['S1', 'S2'],
                                          name='observed_features')
        observed_obs_features = output[1]['observed_features'].view(pd.Series)
        pdt.assert_series_equal(observed_obs_features, expected_obs_features)

        # expected dms, pcoas, and plots returned
        self.assertEqual(set(output[2].keys()), set(['jaccard', 'braycurtis']))
        self.assertEqual(set(output[3].keys()), set(['jaccard', 'braycurtis']))
        self.assertEqual(set(output[4].keys()), set(['jaccard', 'braycurtis']))
        expected_jaccard = skbio.DistanceMatrix([[0, 0.5], [0.5, 0]],
                                                ids=['S1', 'S2'])
        observed_jaccard = output[2]['jaccard'].view(skbio.DistanceMatrix)
        self.assertEqual(observed_jaccard, expected_jaccard)

        # vizard scatter plot returned
        self.assertEqual(output[5].type, Visualization)

    def test_core_metrics_w_replacement(self):
        output = self.core_metrics(table=self.table1,
                                   sampling_depth=2,
                                   metadata=self.metadata,
                                   replacement=True,
                                   n=99)
        # n tables returned
        self.assertEqual(len(output[0]), 99)
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
            observed_table = e.view(pd.DataFrame)
            if observed_table.equals(possible_table1):
                count_possible_table1_observed += 1
            elif observed_table.equals(possible_table2):
                count_possible_table2_observed += 1
            elif observed_table.equals(possible_table3):
                count_possible_table3_observed += 1
            else:
                count_other_table_observed += 1
        # the following assertTrue tests may fail occasionally.
        self.assertTrue(count_possible_table1_observed > 0)
        self.assertTrue(count_possible_table2_observed > 0)
        self.assertTrue(count_possible_table3_observed > 0)
        self.assertEqual(count_other_table_observed, 0)

        # expected alpha vectors returned
        skbio_lt_060_alpha_keys = set(
            ['observed_features', 'pielou_evenness', 'shannon_entropy'])
        skbio_gte_060_alpha_keys = set(
            ['observed_features', 'pielou_e', 'shannon'])
        self.assertTrue(set(output[1].keys()) == skbio_lt_060_alpha_keys or
                        set(output[1].keys()) == skbio_gte_060_alpha_keys)
        observed_obs_features = output[1]['observed_features'].view(pd.Series)
        # note that because n is an odd number, the median for S1 will always
        # be one of the actual values that were observed (as opposed to
        # possibly being 1.5, if n was an even number).
        self.assertTrue(observed_obs_features['S1'] == 1.0 or
                        observed_obs_features['S1'] == 2.0,
                        msg=f"Median value ({observed_obs_features['S1']}) is "
                            "not equal to 1.0 or 2.0.")
        self.assertEqual(observed_obs_features['S2'], 1.0)

        # expected dms, pcoas, and plots returned
        self.assertEqual(set(output[2].keys()), set(['jaccard', 'braycurtis']))
        self.assertEqual(set(output[3].keys()), set(['jaccard', 'braycurtis']))
        self.assertEqual(set(output[4].keys()), set(['jaccard', 'braycurtis']))
        observed_jaccard = output[2]['jaccard'].view(skbio.DistanceMatrix)
        # because n is odd, we should always observe an actual distance
        # between S1 and S2 as the median (as opposed to the mean of two
        # actual values)
        self.assertTrue(observed_jaccard[('S1', 'S2')] in [0.0, 0.5, 1.0],
                        msg=(f"Median value ({observed_jaccard[('S1', 'S2')]})"
                             " is not equal to 0.0, 0.5 or 1.0."))

        # vizard scatter plot returned
        self.assertEqual(output[5].type, Visualization)

    def test_core_metrics_phylogenetic(self):
        output = self.core_metrics(table=self.table1,
                                   phylogeny=self.phylogeny,
                                   sampling_depth=2,
                                   metadata=self.metadata,
                                   replacement=False,
                                   n=10)
        # n tables returned
        self.assertEqual(len(output[0]), 10)
        expected_table = pd.DataFrame(data=[[1.0, 1.0], [0.0, 2.0]],
                                      columns=['F1', 'F2'],
                                      index=['S1', 'S2'])
        for e in output[0].values():
            observed_table = e.view(pd.DataFrame)
            pdt.assert_frame_equal(observed_table, expected_table)

        # expected alpha vectors returned
        skbio_lt_060_alpha_keys = set(
            ['observed_features', 'pielou_evenness', 'shannon_entropy',
             'faith_pd'])
        skbio_gte_060_alpha_keys = set(
            ['observed_features', 'pielou_e', 'shannon', 'faith_pd'])
        self.assertTrue(set(output[1].keys()) == skbio_lt_060_alpha_keys or
                        set(output[1].keys()) == skbio_gte_060_alpha_keys)
        expected_obs_features = pd.Series([2.0, 1.0],
                                          index=['S1', 'S2'],
                                          name='observed_features')
        observed_obs_features = output[1]['observed_features'].view(pd.Series)
        pdt.assert_series_equal(observed_obs_features, expected_obs_features)

        expected_faith_pd = pd.Series([4.0, 3.0],
                                      index=['S1', 'S2'],
                                      name='faith_pd')
        observed_faith_pd = output[1]['faith_pd'].view(pd.Series)
        pdt.assert_series_equal(observed_faith_pd, expected_faith_pd)

        # expected dms, pcoas, and plots returned
        expected_bdiv_keys = set(['jaccard', 'braycurtis',
                                  'unweighted_unifrac', 'weighted_unifrac'])
        self.assertEqual(set(output[2].keys()), expected_bdiv_keys)
        self.assertEqual(set(output[3].keys()), expected_bdiv_keys)
        self.assertEqual(set(output[4].keys()), expected_bdiv_keys)
        expected_jaccard = skbio.DistanceMatrix([[0, 0.5], [0.5, 0]],
                                                ids=['S1', 'S2'])
        observed_jaccard = output[2]['jaccard'].view(skbio.DistanceMatrix)
        self.assertEqual(observed_jaccard, expected_jaccard)

        expected_unweighted_unifrac = skbio.DistanceMatrix(
            [[0, 0.25], [0.25, 0]], ids=['S1', 'S2'])
        observed_unweighted_unifrac = \
            output[2]['unweighted_unifrac'].view(skbio.DistanceMatrix)
        # Floating point error seemingly induced by going from
        # unifrac_binaries=1.4 to unifrac_binaries=1.5 made this necessary
        self.assertEqual(observed_unweighted_unifrac.ids,
                         expected_unweighted_unifrac.ids)
        npt.assert_allclose(observed_unweighted_unifrac.data,
                            expected_unweighted_unifrac.data)

        # vizard scatter plot returned
        self.assertEqual(output[5].type, Visualization)
