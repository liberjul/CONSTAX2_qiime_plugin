# -------------------------------------------------------------------------
# Copyright (c) 2020-2022, JAL&GMNB&GMB
#
# Distributed under the terms of the MIT License.
#
# The full license in the file LICENSE, distributed with this software.
# -------------------------------------------------------------------------
import subprocess
from q2_types.feature_data import (FeatureData, Taxonomy, Sequence, DNAIterator, DNAFASTAFormat)
from qiime2.plugin import Int, Str, Float, Choices, Range
from .plugin_setup import plugin, citations
from ._format_data import _format_ref_db, _detect_format


#python "$CONSTAXPATH"/FormatRefDB.py -d "$DB" -t "$TFILES" -f $FORMAT -p "$CONSTAXPATH"
format = _detect_format(db)
_format_ref_db(db, tfiles, format, dup=False)

print("__________________________________________________________________________\nTraining SINTAX Classifier")

#  "$SINTAXPATH" -makeudb_usearch "${TFILES}/${base}"__UTAX.fasta -output ${TFILES}/sintax.db
cmd = ['vsearch', '-makeudb_usearch', utax_fasta, '-output', sintax_db]
subprocess.run(cmd, check=True)

print("__________________________________________________________________________\nTraining BLAST Classifier")
# makeblastdb -in "${TFILES}/${base}"__RDP_trained.fasta -dbtype nucl -out "${TFILES}/${base}"__BLAST
cmd = ['makeblastdb', '-in', rdp_trained_fasta, '-dbtype', 'nucl', '-out', F'{tf_prefix}__BLAST']
try:
    subprocess.run(cmd, check=True)
except subprocess.CalledProcessError as e:
    print(str(e))
print("__________________________________________________________________________\nTraining RDP Classifier")
#"$RDPPATH" train -o "${TFILES}/." -s "${TFILES}/${base}"__RDP_trained.fasta -t "${TFILES}/${base}"__RDP_taxonomy_trained.txt -Xmx"$MEM"m > rdp_train.out 2>&1
cmd = ['classifier', 'train', '-o', tfiles, '-s', rdp_trained_fasta, '-t', rdp_trained_taxonomy, F'-Xmx{mem}m']
rdp_out = subprocess.run(cmd, capture_output = True)
#Duplicate taxa (a given taxon in more than one higher taxa) occur in some UNITE and SILVA datasets
if "duplicate taxon name" in rdp_out:
    print("RDP training error, redoing with duplicate taxa")
    _format_ref_db(db, tfiles, format, constax_path, dup=True)
    rdp_out = subprocess.run(cmd, capture_output = True)
    if len(rdp.stderr.decode('utf-8')) == 0:
        print("RDP training error overcome, continuing with classification after SINTAX is retrained")
        cmd = ['vsearch', '-makeudb_usearch', utax_fasta, '-output', sintax_db]
        subprocess.run(cmd, check=True)
    else:
        print(rdp.stderr.decode('utf-8'))
        sys.exit(1)

blast_ver = subprocess.run(['blastn', '-version'], capture_output = True).stdout.decode('utf-8').split("\n")[0].split(" ")[1]
with open(F"{tfiles}/training_check.txt", "w") as ofile:
    ofile.write(F"BLAST version {blast_ver}")

plugin.methods.register_function()
