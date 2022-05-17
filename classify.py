# -------------------------------------------------------------------------
# Copyright (c) 2020-2022, JAL&GMNB&GMB

# Distributed under the terms of the MIT License.

# The full license in the file LICENSE, distributed with this software.
# -------------------------------------------------------------------------

import os
import pandas as pd
from q2_types.feature_data import (FeatureData, Taxonomy, Sequence, DNAIterator, DNAFASTAFormat, TSVTaxonomyFormat)
from qiime2.plugin import Int, Str, Float, Choices, Range
from .plugin_setup import plugin, citations
from .train import JSONFormat
from ._format_data import _check_input_names, _split_inputs, _detect_format

def classify(db: DNAFASTAFormat, input: DNAFASTAFormat, training_result: JSONFormat, output_dir: str = 'output', mem: int = 4000, conf: float = 0.8, tax: str = "taxonomy_assignments",
                     nthreads: int = 1, evalue: float = 1., mhits: int = 10, p_iden: float = 0., tf: str = "training_files",
                     isolates: str = "null", iso_qc: int = 75, iso_id: float = 1., hl: str = "null", hl_qc: int = 75, hl_id: float = 1.,
                     conservative: bool = False , consistent: bool = False) -> TSVTaxonomyFormat:
    train_dict = training_result.open()
    format = train_dict["format"]
    db = train_dict['blast_database'].strip("__BLAST")
    if not combine_only:
        formatted_inputs = _check_input_names(input)

        print("__________________________________________________________________________")
        print("Assigning taxonomy to OTU's representative sequences")

        # SINTAX
        sintax_db = F"{tax}/sintax.db"
        cmd = ['vsearch', '-sintax', formatted_inputs, '-db', train_dict['sintax_database'], '-tabbedout',
            F'{tax}/otu_taxonomy.sintax', '-strand', 'both', '-sintax_cutoff', str(conf), '-threads', str(nthreads)]
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
                    cmd = ['blastn', '-query', i, '-db', train_dict['blast_database'], '-num_threads', 'nthreads', '-outfmt', '7', 'qacc', 'sacc', 'evalue', 'bitscore', 'pident', 'qcovs"', '-max_target_seqs', str(mhits)]
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
        cmd = ['classifier', 'classify', '--conf', str(conf), '--format', 'allrank', '--train_propfile', F'{train_dict["rdp_path"]}rRNAClassifier.properties',
            '-o', F'{tax}/otu_taxonomy.rdp', formatted_inputs, F'-Xmx{mem}m']
        subprocess.run(cmd, check=True)

    # Combine classification results
    consensus_Taxonomy = _combine_taxonomy(output_dir, conf, tax, evalue, mhits, p_iden, format, db, tf,
        isolates, iso_qc, iso_id, hl, hl_qc, hl_id, conservative, consistent)
    return consensus_taxonomy

plugin.methods.register_function(
    function=classify,
    inputs={'db' : FeatureData[Sequence],
            'input': FeatureData[Sequence],
            'training_result' : JSONFormat},
    parameters={'mem': Int % Range(1, None),
                'conf' : Float % Range(0., 1., inclusive_end=True, inclusive_start=False),
                'nthreads' : Int % Range(1, None),
                'evalue' : Float % Range(0., None),
                'mhits' : Int % Range(1, None),
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
                'tax': Str},
    outputs=[('classification', FeatureData[Taxonomy])],
    input_descriptions={'db' : 'Database to train classifiers, in FASTA format.',
                        'input': 'Input file in FASTA format containing sequence records to classify.'},
    parameter_descriptions={'mem' : 'Memory available for RDP classification, in MB. Must be in range [1, infinity].',
                            'conf' : 'Classification confidence threshold. Must be in range (0, 1]',
                            'nthreads' : 'Number of threads to use for parallel computing steps. Must be in range [1, infinity].',
                            'evalue' : 'Maximum expect value of BLAST hits to use. Must be in range [0, infinity]',
                            'mhits' : 'Maximum number of BLAST hits to use. Must be in range [1, infinity]',
                            'p_iden' : 'Minimum proportion identity of BLAST hits to use. Must be in range [0, 1]',
                            'tf' : 'Path to which training files will be written',
                            'isolates' : 'FASTA formatted file of isolates to use BLAST against.',
                            'iso_qc': 'Minimum sequence percent query coverage to report isolate matches. Must be in range [0, 100].',
                            'iso_id': 'Minimum aligned sequence proportion identity to report isolate matches. Must be in range [0, 1]',
                            'hl': 'FASTA database file of representative sequences for assignment of high level taxonomy.',
                            'hl_qc': 'Minimum sequence percent query coverage to report high-level taxonomy matches. Must be in range [0, 100].',
                            'hl_id': 'Minimum aligned sequence proportion identity to report high-level taxonomy matches. Must be in range [0, 1]',
                            'conservative': 'Use conservative consensus rule (2 False = False winner).',
                            'consistent': 'Show if the consensus taxonomy is consistent with the real hierarchical taxonomy.',
                            'output_dir': 'Output directory for classifications.',
                            'tax': 'Directory for intermediate taxonomy assignments.'},
    output_descriptions={'classification': 'Taxonomy classifications of query sequences with accompanying statistics and matches to high-level database and/or isolates.'},
    name='CONSTAX2 consensus taxonomy classifier',
    description=(),
    citations=[citations['liber2021constax2']]
    )
