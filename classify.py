# -------------------------------------------------------------------------
# Copyright (c) 2020-2022, JAL&GMNB&GMB

# Distributed under the terms of the MIT License.

# The full license in the file LICENSE, distributed with this software.
# -------------------------------------------------------------------------

import os
import pandas as pd
from q2_types.feature_data import (FeatureData, Taxonomy, Sequence, DNAIterator, DNAFASTAFormat)
from qiime2.plugin import Int, Str, Float, Choices, Range
from .plugin_setup import plugin, citations

from ._format_data import _check_input_names, _split_inputs, _detect_format

def classify(db: DNAFASTAFormat, input: DNAFASTAFormat, output_dir: str = 'output', conf: float = 0.8, tax: str = "taxonomy_assignments",
                     nthreads: int = 1, evalue: float = 1., mhits: int = 10, p_iden: float = 0., tf: str = "training_files",
                     isolates: str = "null", iso_qc: int = 75, iso_id: float = 1., hl: str = "null", hl_qc: int = 75, hl_id: float = 1.,
                     conservative: bool = False , consistent: bool = False) -> pd.DataFrame:
    format = _detect_format(db)
    if ! combine_only:
        formatted_inputs = _check_input_names(input)

        print("__________________________________________________________________________")
        print("Assigning taxonomy to OTU's representative sequences")

        # SINTAX
        sintax_db = F"{tax}/sintax.db"
        cmd = ['vsearch', '-sintax', formatted_inputs, '-db', sintax_db, '-tabbedout',
            F'{tax}/otu_taxonomy.sintax', '-strand', 'both', '-sintax_cutoff', str(conf), '-threads', str(nthreads)"]
        subprocess.run(cmd, check=True)

        cmd = ['sed', '-i', '',  '-e', 's|([0-1][.][0-9]\{2\}|&00|g', F'{tax}/otu_taxonomy.sintax']
        subprocess.run(cmd, check=True)

        # BLAST
        split_queries = _split_inputs(formatted_inputs)
        if isolates != "null":
            formatted_isolates = _check_input_names(isolates, name=F"{tax}/isolates_formatted.fasta")
            iso_blast_db = F'{tax}/{os.path.basename(isolates).split(".")[0]}__BLAST'
            cmd = ['makeblastdb', '-in', F'{tax}/isolates_formatted.fasta', '-dbtype', 'nucl', '-out', iso_blast_db]
            subprocess.run(cmd, check=True)
            iso_buf = ""
        if high_level_db != "null":
            formatted_hl = _check_input_names(high_level_db, name=F"{tax}/hl_formatted.fasta")
            hl_blast_db = F'{tax}/{os.path.basename(high_level_db).split(".")[0]}__BLAST'
            cmd = ['makeblastdb', '-in', F'{tax}/hl_formatted.fasta', '-dbtype', 'nucl', '-out', hl_blast_db]
            subprocess.run(cmd, check=True)
            hl_buf = ""
        with open(F'{tax}/blast.out', 'w') as ofile:
                for i in split_queries:
                    cmd = ['blastn', '-query', i, '-db', F'{tf_prefix}__BLAST', '-num_threads', 'nthreads', '-outfmt', '7', 'qacc', 'sacc', 'evalue', 'bitscore', 'pident', 'qcovs"', '-max_target_seqs', str(mhits)]
                    blast_out = subprocess.run(cmd, capture_output = True).stdout.decode('utf-8')
                    ofile.write(blast_out)
                    if isolates != "null":
                        cmd = ['blastn', '-query', i, '-db', iso_blast_db, '-num_threads', 'nthreads', '-outfmt', '7', 'qacc', 'sacc', 'evalue', 'bitscore', 'pident', 'qcovs"', '-max_target_seqs', '1', '-evalue', '0.00001']
                        iso_buf += subprocess.run(cmd, capture_output = True).stdout.decode('utf-8')
                    if high_level_db != "null":
                        cmd = ['blastn', '-query', i, '-db', hl_blast_db, '-num_threads', 'nthreads', '-outfmt', '7', 'qacc', 'sacc', 'evalue', 'bitscore', 'pident', 'qcovs"', '-max_target_seqs', '1', '-evalue', '0.001']
                        hl_buf += subprocess.run(cmd, capture_output = True).stdout.decode('utf-8')
        if isolates != "null":
            with open(F'{tax}/isolates_blast.out', 'w') as ofile:
                ofile.write(iso_buf)
        if high_level_db != "null":
            with open(F'{tax}/hl_blast.out', 'w') as ofile:
                ofile.write(hl_buf)

        # RDP
        cmd = ['classifier', 'classify', '--conf', str(conf), '--format', 'allrank', '--train_propfile', F'{tf}/rRNAClassifier.properties',
            '-o', F'{tax}/otu_taxonomy.rdp', formatted_inputs, F'-Xmx{mem}m']
        subprocess.run(cmd, check=True)

    # Combine classification results
    consensus_Taxonomy = _combine_taxonomy(output_dir, conf, tax, evalue, mhits, p_iden, format, db, tf,
        isolates, iso_qc, iso_id, hl, hl_qc, hl_id, conservative, consistent)
    return consensus_taxonomy

plugin.methods.register_function(
    function=classify,
    inputs={'db' : FeatureData[Sequence], 'input': FeatureData[Sequence]}
    parameters={'conf' : Float % Range(0., 1., inclusive_end=True),
                'nthreads' : Int,
                'evalue' : Float,
                'mhits' : Int,
                'p_iden' : Float % Range(0., 1., inclusive_end=True),
                'tf' : Str,
                'isolates' : Str,
                'iso_qc': Int % Range(0, 100, inclusive_end=True),
                'iso_id': Float % Range(0., 1., inclusive_end=True),
                'hl': Str,
                'hl_qc': Int % Range(0, 100, inclusive_end=True),
                'hl_id': Float % Range(0., 1., inclusive_end=True),
                'conservative': Bool,
                'consistent': Bool,
                'output_dir': Str,
                'tax': Str}
    outputs=[('classification', FeatureData[Taxonomy])]
    input_descriptions={},
    parameter_descriptions={},
    output_descriptions={},
    name='CONSTAX2 consensus taxonomy classifier',
    description=(),
    citations=[citations['liber2021constax2']]
    )
