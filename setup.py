# ----------------------------------------------------------------------------
# Copyright (c) 2023, Caperaso Lab
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from setuptools import find_packages, setup

import versioneer

setup(
    name="q2-boots",
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    license="BSD-3-Clause",
    packages=find_packages(),
    author="Isaiah Raspet",
    author_email="gregcaperaso@gmail.com",
    description="This is a plugin for bootstraping sample data",
    url="https://github.com/qiime2/q2-boots",
    entry_points={
        "qiime2.plugins": ["q2-boots=q2_boots.plugin_setup:plugin"]
    },
    package_data={
        "q2_boots": ["citations.bib"],
    },
    zip_safe=False,
)
