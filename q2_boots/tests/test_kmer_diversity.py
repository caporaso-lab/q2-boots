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
import pandas.testing as pdt
import skbio


class KmerDiversityTests(TestPluginBase):

    package = 'q2_boots'

    def setUp(self):
        super().setUp()
        self.kmer_diversity = self.plugin.pipelines['kmer_diversity']
        table1 = pd.DataFrame(data=[[1, 1], [0, 4]],
                              columns=['F1', 'F2'],
                              index=['S1', 'S2'])
        self.table1 = qiime2.Artifact.import_data(
            "FeatureTable[Frequency]", table1, view_type=pd.DataFrame)

        sequences = pd.Series(data=[skbio.DNA('GGACCCCTACGCCCATGGTAAACCGACTGGTCGTACGTGA'),  # noqa: E501
                                    skbio.DNA('ACACGGACCTAAGAGCCGACCGCGTACAAAGGCGGGTACGTGCATTGGTTCCGGATCGCCCCGTACATCCGAAGAGCGTC')],  # noqa: E501
                              index=['F1', 'F2'])
        self.sequences1 = qiime2.Artifact.import_data(
            "FeatureData[Sequence]", sequences, view_type=pd.Series)

        metadata = pd.DataFrame(['not', 'of', 'interest'],
                                index=['S1', 'S2', 'S3'],
                                columns=['blank'])
        metadata.index.name = 'sample-id'
        self.metadata = qiime2.Metadata(metadata)

    def test_kmer_diversity_wo_replacement(self):
        output = self.kmer_diversity(table=self.table1,
                                     sequences=self.sequences1,
                                     sampling_depth=2,
                                     metadata=self.metadata,
                                     replacement=False,
                                     n=10)
        # check resampled tables
        self.assertEqual(len(output[0]), 10)
        expected_table = pd.DataFrame(data=[[1.0, 1.0], [0.0, 2.0]],
                                      columns=['F1', 'F2'],
                                      index=['S1', 'S2'])
        for e in output[0].values():
            observed_table = e.view(pd.DataFrame)
            pdt.assert_frame_equal(observed_table, expected_table)

        # check kmer tables
        self.assertEqual(len(output[1]), 10)
        for e in output[1].values():
            observed_table = e.view(pd.DataFrame)
            self.assertTrue(observed_table.shape, (2, 90))

        # expected alpha vectors returned
        skbio_lt_060_alpha_keys = set(
            ['observed_features', 'pielou_evenness', 'shannon_entropy'])
        skbio_gte_060_alpha_keys = set(
            ['observed_features', 'pielou_e', 'shannon'])
        self.assertTrue(set(output[2].keys()) == skbio_lt_060_alpha_keys or
                        set(output[2].keys()) == skbio_gte_060_alpha_keys)
        expected_obs_features = pd.Series([90.0, 65.0],
                                          index=['S1', 'S2'],
                                          name='observed_features')
        observed_obs_features = output[2]['observed_features'].view(pd.Series)
        pdt.assert_series_equal(observed_obs_features, expected_obs_features)

        # expected dms and pcoas returned
        self.assertEqual(set(output[3].keys()), set(['jaccard', 'braycurtis']))
        self.assertEqual(set(output[4].keys()), set(['jaccard', 'braycurtis']))

        # expected values calculated using set operations external to the tests
        expected_jaccard = skbio.DistanceMatrix([[0, 0.27777778],
                                                 [0.27777778, 0]],
                                                ids=['S1', 'S2'])
        observed_jaccard = output[3]['jaccard'].view(skbio.DistanceMatrix)
        pdt.assert_frame_equal(observed_jaccard.to_data_frame(),
                               expected_jaccard.to_data_frame())

        self.assertEqual(output[5].type, Visualization)
