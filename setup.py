# ----------------------------------------------------------------------------
# Copyright (c) 2024, Caporaso Lab (https://cap-lab.bio).
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from setuptools import find_packages, setup

import versioneer

description = ('A QIIME 2 plugin supporting bootstrapped and '
               'rarefaction-based diversity analyses.')

setup(
    name="q2-boots",
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    license="BSD-3-Clause",
    packages=find_packages(),
    author="Isaiah Raspet",
    author_email="caplab@nau.edu",
    description=description,
    url="https://github.com/caporaso-lab/q2-boots",
    entry_points={
        "qiime2.plugins": ["q2-boots=q2_boots.plugin_setup:plugin"]
    },
    package_data={},
    zip_safe=False,
)
