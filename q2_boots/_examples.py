# ----------------------------------------------------------------------------
# Copyright (c) 2024, Caporaso Lab (https://cap-lab.bio).
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import pandas as pd
import qiime2


def table_factory():
    table = pd.DataFrame(data=[[0, 1, 1],
                               [10, 10, 9],
                               [30, 20, 9],
                               [42, 42, 9]],
                         columns=['F1', 'F2', 'F3'],
                         index=['S1', 'S2', 'S3', 'S4'])
    return qiime2.Artifact.import_data(
        "FeatureTable[Frequency]", table, view_type=pd.DataFrame)


def _resample_bootstrap_example(use):
    table = use.init_artifact('table', table_factory)

    resampled_tables, = use.action(
        use.UsageAction(plugin_id='boots',
                        action_id='resample'),
        use.UsageInputs(table=table,
                        sampling_depth=20,
                        n=10,
                        replacement=True),
        use.UsageOutputNames(resampled_tables='bootstrapped_tables')
    )


def _resample_rarefaction_example(use):
    table = use.init_artifact('table', table_factory)

    resampled_tables, = use.action(
        use.UsageAction(plugin_id='boots',
                        action_id='resample'),
        use.UsageInputs(table=table,
                        sampling_depth=20,
                        n=10,
                        replacement=False),
        use.UsageOutputNames(resampled_tables='rarefaction_tables')
    )
