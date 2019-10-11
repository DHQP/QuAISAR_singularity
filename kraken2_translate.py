#!/usr/bin/env python3

#
# Description: Script to convert kraken2 file to label file
#
# Usage: python3 ./kraken2_translate.py -i input_kraken2_file -o output_label_filename
#
# Output location: parameter
#
# Modules required: Biopython must be available in python instance
#
# v1.0 (10/3/2019)
#
# Created by Nick Vlachos (nvx4@cdc.gov)
#

import sys
import glob
import fileinput
import getpass
import argparse
from Bio import Entrez
from time import sleep

# Parse all argument from command line
def parseArgs(args=None):
	parser = argparse.ArgumentParser(description='Script to translate kraken2 file to label file')
	parser.add_argument('-i', '--input', required=True, help='input kraken2 filename')
	parser.add_argument('-o', '--output', required=True, help='output label filename')
	return parser.parse_args()

# Script that will trim fasta files of any sequences that are smaller than the threshold
def get_Taxon_Tree_From_NCBI(taxID):
	Entrez.email = getpass.getuser()
	#Creates the data structure from a pull from entrez nucleotide database using accession id with return type of genbank text mode
	handle = Entrez.efetch(db="taxonomy", id=taxID, mode="text", rettype="xml")
	#Parses the returned output into lines
	result= Entrez.read(handle)
	#Goes through each line until it finds (and prints) the organism name that the accession number represents
	for taxon in result:
		taxid = taxon["TaxId"]
		name = taxon["ScientificName"]
		lineage=["root"]
		print(taxon)
		if "LineageEx" in taxon:
			for t in taxon["LineageEx"]:
				lineage.append(t["ScientificName"])
			lineage.append(name)
		#print("%s\t|\t%s\t|\t%s" % (taxid, name, ";".join(lineage)))
		#print(";".join(lineage))
		return ";".join(lineage)
		#print(' '.join(line.split()[1:]))

# Translate kraken output to a slightly more readable label format.
def translate(input_kraken, output_labels):
	kraken_file=open(input_kraken,'r')
	line=kraken_file.readline().strip()
	tax_tree_dict={}
	mpa_dict={}
	label_lines=[]
	counter=0
	while line != '':
		line_sections = line.split("	")
		contig_id = line_sections[1]
		contig_taxID = line_sections[2]
		if contig_taxID not in tax_tree_dict.keys():
			tax_tree_dict[contig_taxID]=get_Taxon_Tree_From_NCBI(contig_taxID)
		print(str(counter)+":"+contig_id+"	"+tax_tree_dict[contig_taxID])
		label_lines.append(contig_id+"	"+tax_tree_dict[contig_taxID])
		line=kraken_file.readline().strip()
		counter+=1
		sleep(0.005)
	kraken_file.close()
	label_file=open(output_labels, 'w')
	label_file.write("\n".join(label_lines))
	print("Lines:", len(label_lines))
	for line in label_lines:
		print(line)

args = parseArgs()
# Start program
translate(args.input, args.output)
