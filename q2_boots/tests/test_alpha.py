# ----------------------------------------------------------------------------
# Copyright (c) 2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from qiime2.plugin.testing import TestPluginBase


class TestAlphaBootstrap(TestPluginBase):

    package = 'q2_boots'

    def setUp(self):
        super().setUp()
        self.alpha_bootstrap = self.plugin.pipelines['alpha_bootstrap']

    def test_basic(self):
        pass


class TestAlphaPhylogeneticBootstrap(TestPluginBase):

    package = 'q2_boots'

    def setUp(self):
        super().setUp()
        self.alpha_phylogenetic_bootstrap = self.plugin.pipelines[
            'alpha_phylogenetic_bootstrap']

    def test_basic(self):

        pass


class TestAlphaBootstrapRepresentative(TestPluginBase):

    package = 'q2_boots'

    def setUp(self):
        super().setUp()
        self.alpha_bootstrap_representative = self.plugin.pipelines[
            'alpha_bootstrap_representative']

    def test_basic(self):
        pass
