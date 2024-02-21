# ----------------------------------------------------------------------------
# Copyright (c) 2023, QIIME 2 development team.
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


class TestCoreMetrics(TestPluginBase):

    package = 'q2_boots'

    def setUp(self):
        super().setUp()
        self.core_metrics = self.plugin.pipelines['core_metrics']

    def test_basic(self):
        t = Table(data=np.array([[0, 1, 3], [1, 0, 2], [1, 3, 0]]),
                  sample_ids=['S1', 'S2', 'S3'],
                  observation_ids=['O1', 'O2', 'O3'])
        t = Artifact.import_data('FeatureTable[Frequency]', t)
        metadata = pd.DataFrame(['not', 'of', 'interest'],
                                index=['S1', 'S2', 'S3'],
                                columns=['blank'])
        metadata.index.name = 'sample-id'
        metadata = Metadata(metadata)
        output = self.core_metrics(table=t,
                                   sampling_depth=1,
                                   metadata=metadata,
                                   n_jobs=1,
                                   n=10
                                   )
        print(output)
