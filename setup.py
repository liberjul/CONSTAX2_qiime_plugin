# -------------------------------------------------------------------------
# Copyright (c) 2020-2022, JAL&GMNB&GMB
#
# Distributed under the terms of the MIT License.
#
# The full license in the file LICENSE, distributed with this software.
# -------------------------------------------------------------------------

from setuptools import setup, find_packages

import versioneer


setup(
    name="CONSTAX2_qiime_plugin",
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    url="https://qiime2.org",
    license="MIT",
    packages=find_packages(),
    author="Julian A. Liber, Gian M. N. Benucci, and Greg Bonito",
    author_email="julian.liber@duke.edu",
    description="Assign consensus taxonomy to representative sequences. ",
    scripts=['q2_constax/assets/train.py',
             'q2_constax/assets/classify.py'],
    package_data={
        'q2_constax': ['citations.bib'],
        'q2_constax.tests' : ['tests/data/*']
    },
    entry_points={
        "qiime2.plugins":
        ["q2-constax2=q2_constax2.plugin_setup:plugin"]
    },
    zip_safe=False,
)
