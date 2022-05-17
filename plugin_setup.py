# -------------------------------------------------------------------------
# Copyright (c) 2020-2022, JAL&GMNB&GMB
#
# Distributed under the terms of the MIT License.
#
# The full license in the file LICENSE, distributed with this software.
# -------------------------------------------------------------------------

from qiime2.plugin import Plugin, Citations

import CONSTAX2_qiime_plugin

citations = Citations.load('citations.bib', package='q2_CONSTAX2')
plugin = Plugin(
    name='CONSTAX2',
    version=q2_CONSTAX2.__version__,
    website='https://github.com/liberjul/q2_CONSTAX2',
    package='q2_CONSTAX2',
    description=('This QIIME 2 plugin provides consensus'
                'taxonomic classification using SINTAX,'
                'RDP, and BLAST classifiers.'),
    short_description='Plugin for consensus taxonomic classification.',
    citations=[citations['liber2021constax2']]
)
