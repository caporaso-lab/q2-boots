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
import skbio
from io import StringIO


class TestAlphaBootstrap(TestPluginBase):

    package = 'q2_boots'

    def setUp(self):
        super().setUp()
        self.alpha = self.plugin.pipelines['alpha']
        self.alpha_collection = self.plugin.pipelines['alpha_collection']

    def range_check(self, output, div_collection):

        df = None
        i = 0

        all_series = [Artifact.view(x, pd.Series) for x in div_collection]

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

    def test_basic(self):
        t = Table(np.array([[0, 1, 3], [1, 1, 2], [1, 3, 2]]),
                  ['01', '02', '03'],
                  ['S1', 'S2', 'S3'])
        t = Artifact.import_data('FeatureTable[Frequency]', t)
        output = self.alpha(table=t, sampling_depth=1,
                            metric='shannon',
                            random_seed=12,
                            n=10, replacement=True)

        self.assertEqual(len(output), 1)

        output: pd.Series = Artifact.view(output[0], pd.Series)
        self.assertTrue('S1' in output.index)
        self.assertTrue('S2' in output.index)
        self.assertTrue('S3' in output.index)

    def test_range_non_phylo(self):
        t = Table(np.array([[0, 1, 3], [1, 1, 2], [1, 3, 2]]),
                  ['O1', 'O2', 'O3'],
                  ['S1', 'S2', 'S3'])
        t = Artifact.import_data('FeatureTable[Frequency]', t)
        output, = self.alpha(table=t, sampling_depth=1,
                             metric='shannon',
                             random_seed=12,
                             n=10, replacement=True)
        output: pd.Series = Artifact.view(output, pd.Series)

        collection, = self.alpha_collection(
            table=t, sampling_depth=1,
            random_seed=12,
            metric='shannon', n=10, replacement=True
        )

        self.assertTrue(self.range_check(output, collection.values()))

    def test_range_phylo(self):
        with StringIO('(O1:0.3, O2:0.2, O3:0.1, O4:0.2)root;') as f:
            phylogeny = skbio.read(f, format='newick', into=skbio.TreeNode)

        phylogeny = Artifact.import_data(
            'Phylogeny[Rooted]',
            phylogeny
        )
        t = Table(np.array([[0, 1, 3], [1, 1, 2], [1, 3, 2]]),
                  ['O1', 'O2', 'O3'],
                  ['S1', 'S2', 'S3'])
        t = Artifact.import_data('FeatureTable[Frequency]', t)
        output, = self.alpha(table=t, sampling_depth=1,
                             metric='pielou_e',
                             phylogeny=phylogeny,
                             random_seed=12,
                             n=10, replacement=True)
        output: pd.Series = Artifact.view(output, pd.Series)

        collection, = self.alpha_collection(
            table=t, sampling_depth=1,
            phylogeny=phylogeny,
            random_seed=12,
            metric='pielou_e', n=10, replacement=True
        )

        self.assertTrue(self.range_check(output, collection.values()))

    def test_phylo_metric_no_phylo(self):
        t = Table(np.array([[0, 1, 3], [1, 1, 2], [1, 3, 2]]),
                  ['01', '02', '03'],
                  ['S1', 'S2', 'S3'])
        t = Artifact.import_data('FeatureTable[Frequency]', t)

        with self.assertRaisesRegex(ValueError, 'You must use a non-phylogenic metric'):
            self.alpha(table=t, sampling_depth=1,
                       metric='faith_pd',
                       n=10, replacement=True)

    def test_non_phylo_metric_with_phylo(self):
        with StringIO('(O1:0.3, O2:0.2, O3:0.1, O4:0.2)root;') as f:
            phylogeny = skbio.read(f, format='newick', into=skbio.TreeNode)

        phylogeny = Artifact.import_data(
            'Phylogeny[Rooted]',
            phylogeny
        )
        t = Table(np.array([[0, 1, 3], [1, 1, 2], [1, 3, 2]]),
                  ['01', '02', '03'],
                  ['S1', 'S2', 'S3'])
        t = Artifact.import_data('FeatureTable[Frequency]', t)
        # assert no error is thrown, and that phylogeny is just thrown out
        self.alpha(table=t, sampling_depth=1,
                   metric='shannon',
                   random_seed=12,
                   n=10,
                   phylogeny=phylogeny, replacement=True)
        self.assertTrue(True)


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
                                       n=10, replacement=True)

        self.assertEqual(len(output[0]), 10)
        index = ['S1', 'S2', 'S3']

        for series in output[0].values():
            series = Artifact.view(series, pd.Series)
            self.check_samples(index, series)
