# -------------------------------------------------------------------------
# Copyright (c) 2020-2022, JAL&GMNB&GMB
#
# Distributed under the terms of the MIT License.
#
# The full license in the file LICENSE, distributed with this software.
# -------------------------------------------------------------------------

from qiime2.plugin import Plugin, Citations

import CONSTAX2_qiime_plugin

citations = Citations.load('citations.bib', package='CONSTAX2_qiime_plugin')
plugin = Plugin(
    name='CONSTAX2',
    version='v0.0.1', #CONSTAX2_qiime_plugin.__version__,
    website='https://github.com/liberjul/CONSTAX2_qiime_plugin',
    package='CONSTAX2_qiime_plugin',
    description=('This QIIME 2 plugin provides consensus'
                'taxonomic classification using SINTAX,'
                'RDP, and BLAST classifiers.'),
    short_description='Plugin for consensus taxonomic classification.',
    citations=[citations['liber2021constax2']]
)
