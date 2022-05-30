# -------------------------------------------------------------------------
# Copyright (c) 2020-2022, JAL&GMNB&GMB
#
# Includes code originally written by Natalie Vande Pol, November 2017
#
# Distributed under the terms of the MIT License.
#
# The full license in the file LICENSE, distributed with this software.
# -------------------------------------------------------------------------

import sys, os

def _RDP_head_to_UTAX(lineage, format):
	taxa = lineage.split(";")[1:]
	out = ""
	if format != "UNITE" and len(taxa) > 8: # Account for SILVA taxa which have too many ranks to be classified with SINTAX
		taxa = taxa[:8]
	for r in range(len(taxa)):
		if format == "UNITE":
			out = F"{out}{'dkpcofg'[r]}:{taxa[r]},"
		else:
			out = F"{out}{'dkpcofgs'[r]}:{taxa[r]},"

	return out[:-1]
    
def _add_full_lineage(filebase, format):
	print("\n\tAdding Full Lineage\n\n")

	with open(filebase+"__RDP_taxonomy_headers.txt", 'r') as f:
		f1 = f.readlines()
	hash = {} #lineage map

	output_RDP = open(filebase+"__RDP_trained.fasta", 'w')
	output_UTAX = open(filebase+"__UTAX.fasta", 'w')

	for line in f1:
		ID, lineage = line.strip().split("\t")
		ID = ID.strip(">")
		hash[ID] = lineage

	with open(filebase+"__RDP.fasta", 'r') as f:
		f2 = f.readlines()
	for line in f2:
		if line[0] == '>':
			ID = line.strip().replace('>', '')
			try:
				lineage = hash[ID]
			except KeyError:
				print(ID, 'not in taxonomy file')
				sys.exit()
			output_RDP.write(F"{line.strip()}\t{lineage}\n")
			output_UTAX.write(F"{line.strip()};tax={_RDP_head_to_UTAX(lineage, format)};\n")
		else:
			output_RDP.write(line.strip()+"\n")
			output_UTAX.write(line.strip()+"\n")
	output_RDP.close()
	output_UTAX.close()

def _lin_to_tax(file_base, format, dup=False):
	print("\n\tTraining Taxonomy")
	if dup:
		print("\n\tDuplicate taxa being handled with numerical suffices")
	with open(file_base+"__RDP_taxonomy.txt", 'r') as f:
		line = f.readline()
		cols = line.strip().split('\t')[1:] # Split the first line into columns
		hash = {}#taxon name-id map
		ranks = {}#column number-rank map
		hash = {"Root":0}#initiate root rank taxon id map
		for i in range(len(cols)): # Assign ranks based on column headers
			ranks[i] = cols[i]
		root = ['0', 'Root', '-1', '0', 'rootrank']#root rank info
		with open(file_base+"__RDP_taxonomy_trained.txt", 'w') as output_file:
			output_file.write("*".join(root)+"\n")
		ID = 0 #taxon id
		line = f.readline()
		if format == "UNITE":
			name_to_end = {}
			if dup:
				end_name_dict = {}
			while line != "":
				rec_count = 0
				th_buf = "" # taxon header buffer
				output_buf = "" # trained taxonomy buffer
				while line != "" and rec_count < 10000: # Rec count to export when buffer is at 10000 records
					rec_count += 1
					acc = line.split('\t')[0]
					cols = line.strip().split('\t')[1:]
					header = F">{acc}\tRoot"
					for i in range(len(cols)):#iterate each column
						name = []
						for node in cols[:i + 1]:
							if not node == '-':
								name.append(node)
						pName = ";".join(name[:-1])
						depth = len(name)
						name = ";".join(name)
						if name in hash: # Avoid repeated taxonomies
							if name != prev_name:
								prev_name = name
								header += F";{name_to_end[name]}"
							continue
						prev_name = name
						rank = ranks[i]
						if i == 0:
							pName = 'Root'
						pID = hash[pName]#parent taxid
						ID += 1
						hash[name] = ID #add name-id to the map
						end_name = name.split(';')[-1]
						if dup:
							if end_name not in end_name_dict:
								end_name_dict[end_name] = 1
								end_name = F"{end_name}_1"
							else:
								end_name_dict[end_name] += 1
								end_name = F"{end_name}_{end_name_dict[end_name]}"
						header = F"{header};{end_name}"
						name_to_end[name] = end_name
						output_buf = F"{output_buf}{ID}*{end_name}*{pID}*{depth}*{rank}\n"
					th_buf = F"{th_buf}{header}\n"
					line = f.readline()
				with open(file_base+"__RDP_taxonomy_headers.txt", "a+") as taxon_headers:
					taxon_headers.write(th_buf)
				with open(file_base+"__RDP_taxonomy_trained.txt", "a+") as output_file:
					output_file.write(output_buf)
		else:
			name_to_end = {}
			end_name_dict = {}
			while line != "":
				rec_count = 0
				th_buf = ""
				output_buf = ""
				while line != "" and rec_count < 10000:
					rec_count += 1
					acc = line.split('\t')[0]
					cols = line.strip().split('\t')[1:]
					header = F">{acc}\tRoot"
					for i in range(len(cols)):#iterate each column
						name = []
						for node in cols[:i + 1]:
							if not node == '-':
								name.append(node)
						pName = ";".join(name[:-1])
						depth = len(name)
						name = ";".join(name)
						if name in hash:
							if name != prev_name:
								prev_name = name
								header += F";{name_to_end[name]}"
							continue
						prev_name = name
						rank = ranks[i]
						if i == 0:
							pName = 'Root'
						pID = hash[pName]#parent taxid
						ID += 1
						hash[name] = ID #add name-id to the map
						end_name = name.split(';')[-1]
						# Allow for taxa which have more than 1 parent lineage
						if end_name not in end_name_dict:
							end_name_dict[end_name] = 1
							end_name = F"{end_name}_1"
						else:
							end_name_dict[end_name] += 1
							end_name = F"{end_name}_{end_name_dict[end_name]}"
						header = F"{header};{end_name}"
						name_to_end[name] = end_name
						output_buf = F"{output_buf}{ID}*{end_name}*{pID}*{depth}*{rank}\n"
					th_buf = F"{th_buf}{header}\n"
					line = f.readline()
				with open(file_base+"__RDP_taxonomy_headers.txt", "a+") as taxon_headers:
					taxon_headers.write(th_buf)
					print("Headers exported")
				with open(file_base+"__RDP_taxonomy_trained.txt", "a+") as output_file:
					output_file.write(output_buf)
					print("Trained taxonomy exported")

def _convert_lines(line_arr, filter=False, input_name=""):
    if filter and "k__unidentified" in line_arr[0]:
        return ""
    else:
        return F"{unicodedata.normalize('NFKD', line_arr[0]).encode('ASCII', 'ignore').decode().replace(' ', '_')}{check_seq(line_arr[1], input_name, line_arr[0])}"
