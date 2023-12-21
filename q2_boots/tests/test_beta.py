# ----------------------------------------------------------------------------
# Copyright (c) 2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from unittest import TestCase, main
import pandas as pd
from q2_boots.beta import per_cell_average


class TestBeta(TestCase):

    def test_per_cell_median(self):

        a = pd.DataFrame([[1, 1, 1],
                          [1, 1, 1],
                          [1, 1, 1]])
        b = pd.DataFrame([[2, 2, 2],
                          [2, 2, 2],
                          [2, 2, 2]])
        c = pd.DataFrame([[3, 3, 3],
                          [3, 3, 3],
                          [3, 3, 3]])

        result = per_cell_average([a, b, c], 'median')
        exp = pd.DataFrame([[2.0, 2.0, 2.0],
                            [2.0, 2.0, 2.0],
                            [2.0, 2.0, 2.0]])

        pd.testing.assert_frame_equal(exp, result)

    def test_per_cell_mean(self):
        a = pd.DataFrame([[1, 0, 1],
                          [1, 0, 1],
                          [1, 0, 1]])
        b = pd.DataFrame([[2, 2, 2],
                          [2, 2, 2],
                          [2, 2, 2]])
        c = pd.DataFrame([[3, 1, 3],
                          [3, 1, 3],
                          [3, 1, 3]])

        result = per_cell_average([a, b, c], 'mean')
        exp = pd.DataFrame([[2.0, 1.0, 2.0],
                            [2.0, 1.0, 2.0],
                            [2.0, 1.0, 2.0]])

        pd.testing.assert_frame_equal(exp, result)


if __name__ == "__main__":
    main()
