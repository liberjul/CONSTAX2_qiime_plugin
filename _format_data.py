# -------------------------------------------------------------------------
# Copyright (c) 2020-2022, JAL&GMNB&GMB
#
# Includes code originally written by Natalie Vande Pol, November 2017
#
# Distributed under the terms of the MIT License.
#
# The full license in the file LICENSE, distributed with this software.
# -------------------------------------------------------------------------

import sys, os, unicodedata, argparse, glob, time
import numpy as np
from ._lineage import lin_to_tax, add_full_lineage, convert_lines

def _detect_format(db):
    with open(db, "r") as ifile:
    	line = ifile.readline()
    	line_bar_split = line.split("|")
    	if line[0]!=">": # or len(temp0)!= 5:
    		format = "INVALID"
    	elif len(line_bar_split) > 1: # UNITE, because "|" is used in accession
    		if "k__" in line_bar_split[-1]: # kingdom header is defined
    			format="UNITE"
    	elif line.count(";") >= 1: # Silva has ranks divided by ";"
    		format = "SILVA"
    	else:
    		format = "INVALID"
    return format

def _convert_utax_line(utax_taxa):
	r_lets = "dkpcofgs"
	new_taxa = []
	for i in range(min(len(utax_taxa), len(r_lets)) - 1):
		new_taxa.append(F"{r_lets[i]}:{utax_taxa[i]}")
	if utax_taxa[-1].count("_") > 1:
		new_taxa.append(F"{r_lets[-1]}:{utax_taxa[-1]}")
	else:
		new_taxa.append(F"{r_lets[min(len(utax_taxa), len(r_lets)) - 1]}:{utax_taxa[-1]}")
	return ",".join(new_taxa)

def _format_ref_db(db, tf, format, dup=False):
    filename = db
    filename_base = tf + "/" + ".".join(os.path.basename(filename).split(".")[:-1])
    print("\n____________________________________________________________________\nReformatting database\n")
    start = time.process_time()
    fasta = open(filename_base+"__RDP.fasta","w")
    taxon_fn = filename_base+"__RDP_taxonomy.txt"
    taxon = open(taxon_fn,"w")
    print(F"{format} format detected\n")
    if format == "UNITE":
    	taxon.write("Seq_ID\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies\n")

    	num = 0
    	with open(filename) as database:
    		for line in database:
    			if line[0] == ">":
    				#correct umlauts or special letters
    				ascii_line = unicodedata.normalize('NFKD', line).encode('ASCII', 'ignore')
    				temp = ascii_line.decode()[1:].split("|")

    				#RDP files
    				name = str(temp[1])
    				temp2 = temp[4].strip().split("__")
    				to_genus = [ item[:-2] for item in temp2[1:-1] ]

    				if "Incertae_sedis" in to_genus:
    					indices = [i for i,x in enumerate(to_genus) if x == "Incertae_sedis"]
    					for j in indices:
    						if "Incertae_sedis" not in to_genus[j-1]:
    							to_genus[j] = str(to_genus[j-1])+"_Incertae_sedis"
    						else:
    							to_genus[j] = str(to_genus[j-1])
    				if "unidentified" in to_genus:
    					indices = [i for i,x in enumerate(to_genus) if x == "unidentified"]
    					for j in indices:
    						to_genus[j] = "-"

    				if to_genus[0] != "-":
    					species = str(temp2[-1])
    					if "Incertae" in species:
    						species = "unidentified_sp"
    					elif to_genus[-1] not in species:
    						temp=species.split("_")
    						species = temp[0]+"_unidentified_"+temp[1]
    					if species.endswith("sp"):
    						species+= "_"+str(num)
    						num += 1

    					taxonomy = name+"\t"+"\t".join(to_genus)+"\t"+species+"\n"
    					fasta.write(">"+name+"\n")
    					taxon.write(taxonomy)
    					seq = next(database)
    					fasta.write(seq)
    	fasta.close()
    	taxon.close()

    else:
    	max_rank = 1
    	with open(filename, "r") as database:
    		line = database.readline()
    		while line != "":
    			t_list = line.strip().split(";")
    			if len(t_list) > max_rank:
    				max_rank = len(t_list)
    			line = database.readline()
    	taxon.write("Seq_ID\t" + "\t".join([F"Rank_{x}" for x in range(1, max_rank+1)]) + "\n")
    	with open(filename, "r") as database:
    		line = database.readline()
    		while line != "":
    			line = line.replace(" Bacteria;", "?Bacteria;").replace(" Eukaryota;", "?Eukaryota;").replace(" Archaea;", "?Archaea;").replace(" ", "_")
    			ascii_line = unicodedata.normalize('NFKD', line).encode('ASCII', 'ignore')
    			# temp = ascii_line.decode().replace("*", "_").replace("'", "").replace(",", "").replace("Oral_Taxon", "oral_taxon")[1:].split("?")
    			temp = ascii_line.decode().translate(str.maketrans("*,<>", "_   ")).replace("Oral_Taxon", "oral_taxon").replace("'","")[1:].split("?")

    			name = str(temp[0]).split(".")[0]
    			t_list = temp[1].strip().split(";")

    			if "unidentified" in t_list:
    				indices = [i for i,x in enumerate(t_list) if x == "unidentified"]
    				non_unid = [i for i,x in enumerate(t_list) if x != "unidentified"]
    				for j in indices:
    					t_list[j] = "-"
    				if t_list[-1] == "-": # if lowest rank in unidentified, add the lowest identified rank to the lowest taxa
    				    t_list[-1] = t_list[non_unid[-1]] + "_unidentified"

    			if len(t_list) < max_rank:
    				t_list = t_list[:-1] + ["-"]*(max_rank - len(t_list)) + [t_list[-1]] # fill in missing ranks


    			taxonomy = name+"\t"+"\t".join(t_list)+"\n"
    			fasta.write(">"+name+"\n")
    			taxon.write(taxonomy)
    			line = database.readline()
    			seq = ""
    			while line != "" and line[0] != ">":
    				seq += line.strip()
    				line = database.readline()
    			seq = seq.replace("U", "T")
    			fasta.write(seq + "\n")

    	fasta.close()
    	taxon.close()

    print(F"Reference database FASTAs formatted in {time.process_time() - start} seconds...\n")

    os.remove(F"{filename_base}__RDP_taxonomy_trained.txt")
    os.remove(F"{filename_base}__RDP_taxonomy_headers.txt")

    lin_to_tax(filename_base, format, dup)
    add_full_lineage(filename_base, format)

    print("Database formatting complete\n____________________________________________________________________\n\n")

def _check_seq(seq, input_file, otu_name):
    lets = set("ATCGURYSWKMBDHVN\n") #IUPAC nucleotides
    out_seq = seq.upper()
    seq_set = set(out_seq)
    if len(seq_set-lets) > 0:
        raise ValueError(F"Invalid character(s) {seq_set-lets} are present in sequences in input file {input_file}")
    elif len(out_seq) < 16:
        otu_name = otu_name.strip().strip(">")
        warnings.warn(F"Sequence length of {otu_name} is less than 15 nucleotides and may cause this SINTAX error: 'assert failed: m_U.Size == SeqCount'", RuntimeWarning)
    else:
        return out_seq

def _check_input_names(input, name="", filter=False):
    convert_lines_vec = np.vectorize(convert_lines)
    rec_dict={}
    with open(input, "r", encoding='utf-8') as ifile:
        line = ifile.readline()
        while line != "":
            header = line
            line = ifile.readline()
            seq = ""
            while line != "" and line[0] != ">":
                seq = F"{seq}{line}"
                line = ifile.readline()
            rec_dict[header] = seq

    if name == "":
        fname = F"formatted_inputs_{'%06x' % random.randrange(16**6)}.fasta"
        print(fname)
    else:
        fname = name
    rec_array = np.array(list(rec_dict.items()))
    rec_hash = {}
    for i in range(rec_array.shape[0]):
        rec_hash[i] = convert_lines(rec_array[i], filter=filter, input_name=input)

    buffer = "".join(rec_hash.values())
    with open(fname, "w") as ofile:
        ofile.write(buffer)

def _split_inputs(input):
    buffer = ""
    total_rec_count = 0
    file_rec_count = 0
    file_count = 0
    query_files = []
    prefix = input.split(".fasta")[0]
    print("Input FASTA: ", input)
    with open(input, "r") as ifile:
        line = ifile.readline()
        while line != "":
            header = line
            line = ifile.readline()
            seq = ""
            while line != "" and line[0] != ">":
                seq += line.strip()
                line = ifile.readline()
            buffer += F"{header}{seq}\n"
            total_rec_count += 1
            file_rec_count += 1
            if file_rec_count > 99:
                file_rec_count = 0
                with open(F"{prefix}_{file_count:04d}.fasta", "w") as ofile:
                    ofile.write(buffer)
                    query_files.append(F"{prefix}_{file_count:04d}.fasta")
                file_count += 1
                buffer = ""
        if buffer != "":
            with open(F"{prefix}_{file_count:04d}.fasta", "w") as ofile:
                ofile.write(buffer)
                query_files.append(F"{prefix}_{file_count:04d}.fasta")
    return query_files
