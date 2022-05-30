# -------------------------------------------------------------------------
# Copyright (c) 2020-2022, JAL&GMNB&GMB
#
# Distributed under the terms of the MIT License.
#
# The full license in the file LICENSE, distributed with this software.
# -------------------------------------------------------------------------

import os
from q2_types.feature_data import (FeatureData, Taxonomy, Sequence, DNAIterator, DNAFASTAFormat, TSVTaxonomyFormat)

from q2_constax.assets import train
from . import CONSTAXTests

class TrainTests(CONSTAXTests):

    def setUp(self):
        super().setUp()
        self.db_path_unite = self.get_data_path('unite_test_ref.fasta')
        database_unite = Artifact.import_data('FeatureData[Sequence]', self.db_path_unite)



assert os.path.isfile(self.trained_dict_unite['sintax_database'])
assert os.path.isfile(self.trained_dict_silva['sintax_database'])
