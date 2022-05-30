# -------------------------------------------------------------------------
# Copyright (c) 2020-2022, JAL&GMNB&GMB

# Distributed under the terms of the MIT License.

# The full license in the file LICENSE, distributed with this software.
# -------------------------------------------------------------------------

def fasta_select_by_keyword(input, output, keyword):
    rec_dict = {}
    with open(input, "r") as ifile:
        line = ifile.readline()
        while line != "":
            header = line
            line = ifile.readline()
            seq = ""
            while line != "" and line[0] != ">":
                seq += line.strip()
                line = ifile.readline()
            rec_dict[header] = seq
    with open(output, "w") as ofile:
        for rec in rec_dict.keys():
            if keyword in rec:
                ofile.write(F"{rec}{rec_dict[rec]}\n")
