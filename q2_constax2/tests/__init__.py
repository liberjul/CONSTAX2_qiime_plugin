# -------------------------------------------------------------------------
# Copyright (c) 2020-2022, JAL&GMNB&GMB
#
# Distributed under the terms of the MIT License.
#
# The full license in the file LICENSE, distributed with this software.
# -------------------------------------------------------------------------

import tempfile

from qiime2.plugin.testing import TestPluginBase

class CONSTAXTests(TestPluginBase):
    def setUp(self):
        try:
            from q2_constax2.plugin_setup import plugin
        except ImportError:
            self.fail("Could not import plugin object.")

        self.plugin = plugin

        self.temp_dir = tempfile.TemporaryDirectory(
            prefix='q2-constax2-test-temp-')

        filterwarnings('ignore', UserWarning)
