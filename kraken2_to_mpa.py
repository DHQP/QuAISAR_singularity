#!/usr/bin/env python3

#
# Description: Script to convert kraken2 file to mpa file
#
# Usage: python3 ./kraken2_report_from_kraken.py -i input_kraken2_file -o output_mpa_filename
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

# Parse all arguments from command line
def parseArgs(args=None):
	parser = argparse.ArgumentParser(description='Script to turn kraken2 file to mpa')
	parser.add_argument('-i', '--input', required=True, help='input kraken2 filename')
	parser.add_argument('-o', '--output', required=True, help='output mpa filename')
	return parser.parse_args()

# Script that will retrieve the taxonomic lineage of a sample from entrez
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
		for t in taxon["LineageEx"]:
			lineage.append(t["ScientificName"])
		lineage.append(name)
		#print("%s\t|\t%s\t|\t%s" % (taxid, name, ";".join(lineage)))
		return ";".join(lineage)
		#print(' '.join(line.split()[1:]))

# Translate kraken output to a slightly more readable label format.
def translate(input_kraken, output_labels):
	kraken=open(input_kraken,'r')
	line=kraken.readline().strip()
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
		line=kraken.readline().strip()
		counter+=1
	kraken.close()
	label_file=open(output_labels, 'w')
	label_file.write("\n".join(label_lines))
	print("Lines:", len(label_lines))
	for line in label_lines:
		print(line)

#translate(sys.argv[1], sys.argv[2])

# Script that will retrieve taxonomic info from entrez and convert it into mpa format
def get_mpa_string_From_NCBI(taxID):
	Entrez.email = getpass.getuser()
	#Creates the data structure from a pull from entrez nucleotide database using accession id with return type of genbank text mode
	handle = Entrez.efetch(db="taxonomy", id=taxID, mode="text", rettype="xml")
	#Parses the returned output into lines
	result= Entrez.read(handle)
	#Will have to change this region to be able to handle more descriptive lists downstream
	recognized_ranks={"superkingdom":"d", "kingdom":"k", "phylum":"p", "class":"c", "order":"o", "family":"f", "genus":"g", "species":"s", "species_group":"x"}
	for entry in result:
		#taxid = entry["Rank"]
		#print(entry)
		mpa_string=""
		if "LineageEx" in entry:
			for r in entry["LineageEx"]:
				#print(r)
				#print(r["Rank"])
				if r["Rank"] in recognized_ranks.keys():
					current_rank=recognized_ranks[r["Rank"]]
				#else:
					#current_rank="-"
					current_taxa=(r["ScientificName"])
					rank_and_taxa=current_rank+"__"+current_taxa
					#print(rank_and_taxa)
					mpa_string+=rank_and_taxa+"|"
		if entry["Rank"] in recognized_ranks.keys() or entry["Rank"] == "no rank":
			if entry["Rank"] == "no rank":
				current_rank="-"
			else:
				current_rank=recognized_ranks[entry["Rank"]]
		#else:
			#current_rank="-"
			current_taxa=(entry["ScientificName"])
			rank_and_taxa=current_rank+"__"+current_taxa
			#print(rank_and_taxa)
			mpa_string+=rank_and_taxa
		#else:
		#	print(entry["Rank"])
		return(mpa_string)

# Sorts mpas to put similar entries together, easier for humans to look at
def organize_mpas(input_kraken, output_mpa):
	kraken=open(input_kraken,'r')
	line=kraken.readline().strip()
	mpa_dict={}
	mpa_counts={}
	counter=0
	while line != '':
		line_sections = line.split("	")
		contig_id = line_sections[1]
		contig_taxID = line_sections[2]
		if contig_taxID not in mpa_dict.keys():
			print("Adding", contig_taxID)
			mpa_dict[contig_taxID]=get_mpa_string_From_NCBI(contig_taxID)
			mpa_counts[contig_taxID]=1
		else:
			if contig_taxID not in mpa_counts.keys():
				#print("Doesnt exist???", contig_taxID)
				mpa_counts[contig_taxID]=1
			else:
				#print("Incrementing:", contig_taxID)
				mpa_counts[contig_taxID]+=1
		line=kraken.readline().strip()
		sleep(0.005)
	kraken.close()
	mpa_taxon_counts={}
	#print("mpa_dict length:", len(mpa_dict))
	for key in mpa_dict.keys():
		#print(key, mpa_dict[key])
		taxons=mpa_dict[key].split("|")
		if taxons is not None:
			#print(taxons)
			for i in range(0, len(taxons)):
				#print(taxons[0:i+1])
				if taxons[i] != "" and taxons[i][0:1] != "-":
					#print(taxons[i])
					if "|".join(taxons[0:i+1]) in mpa_taxon_counts:
						#print("Incrementing", "|".join(taxons[0:i+1]), "from", mpa_taxon_counts["|".join(taxons[0:i+1])], "to",  mpa_taxon_counts["|".join(taxons[0:i+1])]+mpa_counts[key])
						mpa_taxon_counts["|".join(taxons[0:i+1])]+=mpa_counts[key]
					else:
						#print("Creating", "|".join(taxons[0:i+1]), "at 1")
						mpa_taxon_counts["|".join(taxons[0:i+1])]=mpa_counts[key]
		else:
			print("Taxons is none")

	print("mpa_taxon_counts length:", len(mpa_taxon_counts))
	for key in sorted(mpa_taxon_counts.keys()):
		print(key, mpa_taxon_counts[key])






	#print("mpa_counts length:", len(mpa_counts))
	#for key in mpa_counts.keys():
	#	print(key, mpa_counts[key])


#get_mpa_string_From_NCBI(470, blank_dick)

args = parseArgs()
organize_mpas(args.input, args.output)
