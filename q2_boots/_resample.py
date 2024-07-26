# ----------------------------------------------------------------------------
# Copyright (c) 2024, Caporaso Lab (https://cap-lab.bio).
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

def resample(ctx, table, sampling_depth, n, replacement):
    rarefy_action = ctx.get_action('feature_table', 'rarefy')
    resampled_tables = []

    for i in range(n):
        resampled_table = rarefy_action(table=table,
                                        sampling_depth=sampling_depth,
                                        with_replacement=replacement)[0]
        resampled_tables.append(resampled_table)

    return {f'resampled-table-{i}': t for i, t in enumerate(resampled_tables)}
