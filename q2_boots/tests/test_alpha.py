# ----------------------------------------------------------------------------
# Copyright (c) 2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from qiime2.plugin.testing import TestPluginBase
from biom.table import Table
import numpy as np
import pandas as pd
from qiime2 import Artifact


class TestAlphaBootstrap(TestPluginBase):

    package = 'q2_boots'

    def setUp(self):
        super().setUp()
        self.alpha_bootstrap = self.plugin.pipelines['alpha']

    def test_basic(self):
        t = Table(np.array([[0, 1, 3], [1, 1, 2], [1, 3, 2]]),
                  ['01', '02', '03'],
                  ['S1', 'S2', 'S3'])
        t = Artifact.import_data('FeatureTable[Frequency]', t)
        output = self.alpha_bootstrap(table=t, sampling_depth=1,
                                      metric='shannon',
                                      n=10)

        self.assertEqual(len(output), 1)

        output: pd.Series = Artifact.view(output[0], pd.Series)
        self.assertTrue('S1' in output.index)
        self.assertTrue('S2' in output.index)
        self.assertTrue('S3' in output.index)


class TestAlphaPhylogeneticBootstrap(TestPluginBase):

    package = 'q2_boots'

    def setUp(self):
        super().setUp()
        self.alpha_phylogenetic_bootstrap = self.plugin.pipelines[
            'alpha_collection']

    def test_basic(self):

        pass


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
        t = Artifact.import_data('FeatureTable[Frequency]', t)
        output = self.alpha_collection(table=t, sampling_depth=1,
                                       metric='shannon',
                                       n=10)

        self.assertEqual(len(output[0]), 10)
        index = ['S1', 'S2', 'S3']

        for series in output[0].values():
            series = Artifact.view(series, pd.Series)
            self.check_samples(index, series)
