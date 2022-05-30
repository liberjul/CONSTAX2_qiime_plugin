# -------------------------------------------------------------------------
# Copyright (c) 2020-2022, JAL&GMNB&GMB
#
# Distributed under the terms of the MIT License.
#
# The full license in the file LICENSE, distributed with this software.
# -------------------------------------------------------------------------
import subprocess, os, json

from q2_types.feature_data import (FeatureData, Taxonomy, Sequence, DNAIterator, DNAFASTAFormat)
from qiime2.plugin import Int, Str, Float, Choices, Range, TextFileFormat, SemanticType, BinaryFileFormat
import qiime2.plugin.model as model
from .plugin_setup import plugin, citations
from ._format_data import _format_ref_db, _detect_format

# https://github.com/qiime2/q2-types/issues/49
# From https://github.com/qiime2/q2-feature-classifier/blob/master/q2_feature_classifier/_taxonomic_classifier.py

CONSTAXTaxonomicClassifier = SemanticType('CONSTAXTaxonomicClassifier')

class JSONFormat(model.TextFileFormat):
    def sniff(self):
        with self.open() as fh:
            try:
                json.load(fh)
                return True
            except json.JSONDecodeError:
                pass
        return False

class CONSTAXTaxonomicClassifierDirFmt(model.DirectoryFormat):
    training_result = model.File('training_result.json', format=JSONFormat)


@plugin.register_transformer
def _1(fmt: JSONFormat) -> dict:
    with fmt.open() as fh:
        return json.load(fh)


@plugin.register_transformer
def _2(fmt: CONSTAXTaxonomicClassifierDirFmt) -> dict:
    with fmt.training_result.open() as fh:
        return json.load(fh)

@plugin.register_transformer
def _3(data: dict) -> JSONFormat:
    result = JSONFormat()
    with result.open() as fh:
        json.dump(data, fh)
    return result

def train(db : DNAFASTAFormat, tf : str, mem : int) -> dict:
    db.file.view(DNAFASTAFormat)
    training_dict = {'sintax_database' : F'{tf}/sintax.db',
                     'blast_database' : F'{tf}/{db_base}__BLAST',
                     'rdp_path' : F'{tf}/'}
    #python "$CONSTAXPATH"/FormatRefDB.py -d "$DB" -t "$TFILES" -f $FORMAT -p "$CONSTAXPATH"
    format = _detect_format(db)
    training_dict['format' : format]
    _format_ref_db(db, tf, format, dup=False)
    db_base = os.path.basename(db).split(".")[0]
    print("__________________________________________________________________________\nTraining SINTAX Classifier")

    #  "$SINTAXPATH" -makeudb_usearch "${TFILES}/${base}"__UTAX.fasta -output ${TFILES}/sintax.db
    cmd = ['vsearch', '-makeudb_usearch', F'{tf}/{db_base}__UTAX.fasta', '-output', F'{tf}/sintax.db']
    subprocess.run(cmd, check=True)

    print("__________________________________________________________________________\nTraining BLAST Classifier")
    # makeblastdb -in "${TFILES}/${base}"__RDP_trained.fasta -dbtype nucl -out "${TFILES}/${base}"__BLAST
    cmd = ['makeblastdb', '-in', F'{tf}/{db_base}__RDP_trained.fasta', '-dbtype', 'nucl', '-out', F'{tf}/{db_base}__BLAST']
    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as e:
        print(str(e))
    print("__________________________________________________________________________\nTraining RDP Classifier")
    #"$RDPPATH" train -o "${TFILES}/." -s "${TFILES}/${base}"__RDP_trained.fasta -t "${TFILES}/${base}"__RDP_taxonomy_trained.txt -Xmx"$MEM"m > rdp_train.out 2>&1
    cmd = ['classifier', 'train', '-o', F'{tf}/.', '-s', F'{tf}/{db_base}__RDP_trained.fasta', '-t', F'{tf}/{db_base}__RDP_taxonomy_trained.txt', F'-Xmx{mem}m']
    rdp_out = subprocess.run(cmd, capture_output = True)
    #Duplicate taxa (a given taxon in more than one higher taxa) occur in some UNITE and SILVA datasets
    if "duplicate taxon name" in rdp_out:
        print("RDP training error, redoing with duplicate taxa")
        _format_ref_db(db, tfiles, format, constax_path, dup=True)
        rdp_out = subprocess.run(cmd, capture_output = True)
        if len(rdp.stderr.decode('utf-8')) == 0:
            print("RDP training error overcome, continuing with classification after SINTAX is retrained")
            cmd = ['vsearch', '-makeudb_usearch', F'{tf}/{db_base}__UTAX.fasta', '-output', F'{tf}/sintax.db']
            subprocess.run(cmd, check=True)
        else:
            raise RuntimeError("RDP training with duplicate taxa failed:\n" + rdp.stderr.decode('utf-8'))

    blast_ver = subprocess.run(['blastn', '-version'], capture_output = True).stdout.decode('utf-8').split("\n")[0].split(" ")[1]
    training_dict['blast_version'] = blast_ver

    rdp_version = subprocess.run(['conda', 'list', 'rdptools'], capture_output = True).stdout.decode('utf-8').split("\n")[-2].split()[1]
    with open(F"{tfiles}/rRNAClassifier.properties", "w") as ofile:
        ofile.write(F"# Sample ResourceBundle properties file\nbergeyTree=bergeyTrainingTree.xml\n\nprobabilityList=genus_wordConditionalProbList.txt\n\nprobabilityIndex=wordConditionalProbIndexArr.txt\n\nwordPrior=logWordPrior.txt\n\nclassifierVersion=RDP Naive Bayesian rRNA Classifier Version {rdp_version}")

    json.dump(training_dict, F"{tfiles}/training_result.json")
    training_result = JSONFormat(F"{tfiles}/training_result.json")
    return training_result

plugin.register_semantic_types(CONSTAXTaxonomicClassifier)
plugin.register_semantic_type_to_format(CONSTAXTaxonomicClassifier,
    artifact_format = CONSTAXTaxonomicClassifierDirFmt)

plugin.methods.register_function(
    function=train,
    inputs={'db' : FeatureData[Sequence]},
    parameters={'mem' : Int % Range(0, None),
                'tf' : Str},
    outputs=[('training_result', CONSTAXTaxonomicClassifier)],
    input_descriptions={'db' : 'Database to train classifiers, in FASTA format.'},
    parameter_descriptions={'mem' : 'Memory available for RDP classification, in MB. Must be in range [1, infinity].',
                            'tf' : 'Path to which training files will be written'},
    output_descriptions={'training_result': 'JSON with attributes describing model training to allow for classification.'},
    name='CONSTAX2 consensus taxonomy classifier',
    description='Function to train the CONSTAX classifiers on a reference database',
    citations=[citations['liber2021constax2']]
    )
