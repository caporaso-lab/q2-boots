# ----------------------------------------------------------------------------
# Copyright (c) 2024, Caporaso Lab (https://cap-lab.bio).
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import pandas as pd

import qiime2
from qiime2.plugin.testing import TestPluginBase


class ResampleTests(TestPluginBase):
    package = 'q2_boots.tests'

    def setUp(self):
        super().setUp()
        self.resample_pipeline = self.plugin.pipelines['resample']

        table1 = pd.DataFrame(data=[[0, 1], [1, 1], [3, 2]],
                              columns=['F1', 'F2'],
                              index=['S1', 'S2', 'S3'])
        self.table_artifact1 = qiime2.Artifact.import_data(
            "FeatureTable[Frequency]", table1, view_type=pd.DataFrame
        )

        table2 = pd.DataFrame(data=[[0, 1, 1],
                                    [10, 10, 9],
                                    [30, 20, 9],
                                    [42, 42, 9]],
                              columns=['F1', 'F2', 'F3'],
                              index=['S1', 'S2', 'S3', 'S4'])
        self.table_artifact2 = qiime2.Artifact.import_data(
            "FeatureTable[Frequency]", table2, view_type=pd.DataFrame
        )

        table3 = pd.DataFrame(data=[[1, 1]],
                              columns=['F1', 'F2'],
                              index=['S1'])
        self.table_artifact3 = qiime2.Artifact.import_data(
            "FeatureTable[Frequency]", table3, view_type=pd.DataFrame
        )

    def test_resample_w_replacement_filters_samples(self):
        self._resample_filters_sample(replacement=True)

    def test_resample_wo_replacement_filters_samples(self):
        self._resample_filters_sample(replacement=False)

    def test_expected_sampling_depth_w_replacement(self):
        self._expected_sampling_depth(replacement=True)

    def test_expected_sampling_depth_wo_replacement(self):
        self._expected_sampling_depth(replacement=False)

    def test_expected_n_tables(self):
        obs_tables, = self.resample_pipeline(table=self.table_artifact1,
                                             sampling_depth=1,
                                             n=4,
                                             replacement=True)
        self.assertEqual(len(obs_tables), 4)

        obs_tables, = self.resample_pipeline(table=self.table_artifact1,
                                             sampling_depth=1,
                                             n=2,
                                             replacement=True)
        self.assertEqual(len(obs_tables), 2)

    def test_w_replacement(self):
        obs_tables, = self.resample_pipeline(table=self.table_artifact3,
                                             sampling_depth=2,
                                             n=100,
                                             replacement=True)
        # if sampling with replacement from a sample with 2 unique features
        # that have one observation each, we should observe a resampled
        # table with only one feature 50% of the time. the probability of not
        # seeing a table with only one feature in 100 resample tables is 1e-30,
        # so intermittent failure of this test is possible but should be
        # extremely infrequent
        fewer_than_two_unique_features_ever_observed = False
        for obs_table in obs_tables.values():
            obs_table = obs_table.view(pd.DataFrame)
            if len(obs_table.columns) < 2:
                fewer_than_two_unique_features_ever_observed = True
        self.assertTrue(fewer_than_two_unique_features_ever_observed)

    def test_wo_replacement(self):
        obs_tables, = self.resample_pipeline(table=self.table_artifact3,
                                             sampling_depth=2,
                                             n=100,
                                             replacement=False)
        # if sampling with replacement from a sample with 2 unique features
        # that have one observation each, we should never observe a resampled
        # feature table with only one feature. confirm that over 100
        # iterations we always have two features in the resampled table.
        exactly_two_features_always_observed = True
        for obs_table in obs_tables.values():
            obs_table = obs_table.view(pd.DataFrame)
            if len(obs_table.columns) != 2:
                exactly_two_features_always_observed = False
        self.assertTrue(exactly_two_features_always_observed)

    # test helper functions
    def _expected_sampling_depth(self, replacement):
        obs_tables, = self.resample_pipeline(table=self.table_artifact2,
                                             sampling_depth=1,
                                             n=10,
                                             replacement=replacement)
        for obs_table in obs_tables.values():
            obs_table = obs_table.view(pd.DataFrame)
            self.assertEqual(list(obs_table.sum(axis=1)),
                             [1., 1., 1., 1.])

        obs_tables, = self.resample_pipeline(table=self.table_artifact2,
                                             sampling_depth=2,
                                             n=10,
                                             replacement=replacement)
        for obs_table in obs_tables.values():
            obs_table = obs_table.view(pd.DataFrame)
            self.assertEqual(list(obs_table.sum(axis=1)),
                             [2., 2., 2., 2.])

        obs_tables, = self.resample_pipeline(table=self.table_artifact2,
                                             sampling_depth=50,
                                             n=10,
                                             replacement=replacement)
        for obs_table in obs_tables.values():
            obs_table = obs_table.view(pd.DataFrame)
            self.assertEqual(list(obs_table.sum(axis=1)),
                             [50., 50.])

    def _resample_filters_sample(self, replacement):
        with self.assertRaisesRegex(ValueError, "no samples or features"):
            _ = self.resample_pipeline(table=self.table_artifact1,
                                       sampling_depth=6,
                                       n=10,
                                       replacement=replacement)

        obs_tables, = self.resample_pipeline(table=self.table_artifact1,
                                             sampling_depth=5,
                                             n=10,
                                             replacement=replacement)
        for obs_table in obs_tables.values():
            obs_table = obs_table.view(pd.DataFrame)
            sids = list(obs_table.index)
            self.assertEqual(sids, ['S3'])

        obs_tables, = self.resample_pipeline(table=self.table_artifact1,
                                             sampling_depth=2,
                                             n=10,
                                             replacement=replacement)
        for obs_table in obs_tables.values():
            obs_table = obs_table.view(pd.DataFrame)
            sids = list(obs_table.index)
            self.assertEqual(sids, ['S2', 'S3'])

        obs_tables, = self.resample_pipeline(table=self.table_artifact1,
                                             sampling_depth=1,
                                             n=10,
                                             replacement=replacement)
        for obs_table in obs_tables.values():
            obs_table = obs_table.view(pd.DataFrame)
            sids = list(obs_table.index)
            self.assertEqual(sids, ['S1', 'S2', 'S3'])
