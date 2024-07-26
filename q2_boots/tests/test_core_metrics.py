# ----------------------------------------------------------------------------
# Copyright (c) 2024, Caporaso Lab (https://cap-lab.bio).
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from qiime2.plugin.testing import TestPluginBase
from qiime2.metadata import Metadata
from biom.table import Table
from qiime2 import Artifact
import numpy as np
import pandas as pd
import skbio
from io import StringIO


class TestCoreMetrics(TestPluginBase):

    package = 'q2_boots'

    def setUp(self):
        super().setUp()
        self.core_metrics = self.plugin.pipelines['core_metrics']

    def test_basic(self):
        t = Table(data=np.array([[1000, 1920, 3451], [4536, 1552, 6521],
                                 [1634, 1634, 6721]]),
                  sample_ids=['S1', 'S2', 'S3'],
                  observation_ids=['O1', 'O2', 'O3'])
        t = Artifact.import_data('FeatureTable[Frequency]', t)
        metadata = pd.DataFrame(['not', 'of', 'interest'],
                                index=['S1', 'S2', 'S3'],
                                columns=['blank'])
        metadata.index.name = 'sample-id'
        metadata = Metadata(metadata)
        output = self.core_metrics(table=t,
                                   sampling_depth=500,
                                   metadata=metadata,
                                   n_jobs=1,
                                   replacement=True,
                                   n=10
                                   )
        self.assertEqual(len(output[0]), 10)

        for table in output[0].values():
            table = Artifact.view(table, pd.DataFrame)
            self.assertEqual(table.shape, (3, 3))

    def test_phylogeny(self):
        with StringIO('(O1:0.3, O2:0.2, O3:0.1, O4:0.2)root;') as f:
            phylogeny = skbio.read(f, format='newick', into=skbio.TreeNode)

        phylogeny = Artifact.import_data(
            'Phylogeny[Rooted]',
            phylogeny
        )
        metadata = pd.DataFrame(['not', 'of', 'interest'],
                                index=['S1', 'S2', 'S3'],
                                columns=['blank'])
        metadata.index.name = 'sample-id'
        metadata = Metadata(metadata)
        t = Table(data=np.array([[1000, 1920, 3451], [2536, 1552, 1521],
                                 [1634, 1634, 6721]]),
                  sample_ids=['S1', 'S2', 'S3'],
                  observation_ids=['O1', 'O2', 'O3'])
        t = Artifact.import_data('FeatureTable[Frequency]', t)
        output = self.core_metrics(table=t,
                                   sampling_depth=500,
                                   metadata=metadata,
                                   n_jobs=1,
                                   n=10,
                                   phylogeny=phylogeny,
                                   replacement=True,
                                   )
        self.assertEqual(len(output[0]), 10)

        for table in output[0].values():
            table = Artifact.view(table, pd.DataFrame)
            self.assertEqual(table.shape, (3, 3))
