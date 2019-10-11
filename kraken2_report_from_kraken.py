#!/usr/bin/env python3

#
# Description: Script to convert kraken2 file to list file
#
# Usage: python3 ./kraken2_report_from_kraken.py -i input_kraken2_file -o output_list_filename
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
from Bio import Entrez
import argparse
from operator import attrgetter

# Parse all arguments from command line
def parseArgs(args=None):
	parser = argparse.ArgumentParser(description='Script to convert kraken2 file to list')
	parser.add_argument('-i', '--input', required=True, help='input kraken2 filename')
	parser.add_argument('-o', '--output', required=True, help='output list filename')
	return parser.parse_args()

# Class to represent a taxonomic node
class taxon_Node:
	#initializing the variables
	name = ""
	count = 0
	ID = 0
	parent = None
	level = "-"
	children = []

	#defining constructor
	def __init__(self, taxonName, taxonCount, taxonParent, taxonChildren, taxonID, taxonLevel):
		self.name = taxonName
		self.count = taxonCount
		self.parent = taxonParent
		if taxonChildren is not None:
			self.children = taxonChildren
		self.ID = taxonID
		self.level = taxonLevel

	def __eq__(self, other):
		return self.count == other.count

	def __lt__(self, other):
		return self.count < other.count

	#defining class methods
	def showName(self):
		#print("A")
		print("Name:", self.name)
		#print("B")

	def getName(self):
		return self.name

	def showCount(self):
		print("Counts:", self.count)

	def getCount(self):
		return self.count

	def showParentName(self):
		if self.parent is not None:
			print("Parent:")
			self.parent.showName()
		else:
			print("No parent assigned")

	def getParent(self):
		return self.parent

	def showtaxID(self):
		print("taxID:", self.ID)

	def getTaxID(self):
		return taxID

	def showChildren(self):
		if self.children is None or len(self.children) == 0:
			print("No Children")
		else:
			print("Children:", len(self.children))
			for child in self.children:
				#print("Child:")
				child.showName()
				#print("dlihc")

	def getChildren(self):
		return self.children

	def addCounts(self, newReads):
		self.count += newReads

	def addChild(self, newChild):
		#print("Adding")
		#newChild.showName()
		#newChild.showtaxID()
		if self.getChildCount() == 0:
			#newChild.setParent(self)
			self.children=[newChild]
			#print("Children:", len(self.children), self.children[0].getName())
		else:
			#newChild.setParent(self)
			self.children.append(newChild)
			#print("Children:", len(self.children), self.children[0].getName(), self.children[len(self.children)-1].getName())
		#print("End Adding")

	def isChild_name(self, checkName):
		if self.children is not None:
			for child in self.children:
				if child.getName() == checkName:
					return True
		return False

	def isChild_taxon(self, checkName):
		if self.children is not None:
			for child in self.children:
				if child.getTaxID() == checkName:
					return True
		return False

	def getUnpassableCounts(self):
		if self.getChildCount() == 0:
			return self.count
		else:
			unpassable=self.count
			for child in self.children:
				unpassable-=child.getCount()
			return unpassable

	def find_taxID(self, checkName):
		if self.getChildCount() > 1:
			for child in self.children:
				print(len(self.children), child.getTaxID(), checkName)
				if child.getTaxID() == checkName:
					return child
				else:
					child.find(checkName)
		return None

	def find_name(self, checkName):
		if self.name == checkName:
			print("Found")
			return self
		elif self.getChildCount() > 0:
			for child in self.children:
				print(len(self.children), child.getName(), checkName)
				temp_result = child.find_name(checkName)
				if temp_result is not None:
					return temp_result

		else:
			return None

	def setParent(self, newParent):
		self.parent = newParent

	def getChildCount(self):
		if self.children is None or len(self.children) == 0:
			return 0
		else:
			return len(self.children)

	def showLevel(self):
		print("level:", self.level)


	def print(self):
		print("Name:", self.name +"\nCounts:", str(self.count) +"\nID:", str(self.ID))
		print(self.getChildCount())
		self.showParent()
		self.showLevel()
		#print("End Print")

	def print_All(self):
		print("Name:", self.name +"\nCounts:", str(self.count) +"\nID:", str(self.ID))
		print("#Children:", self.getChildCount())
		self.showParentName()
		self.showLevel()
		if self.getChildCount() is not 0:
			self.sort_Children()
			for child in self.children:
				child.print_All()

	def sort_Children(self):
		if self.getChildCount() > 1:
			print("Testing rearrangement")
			self.showChildren()
			self.children.sort(key=lambda x: x.count, reverse = True)
			self.showChildren()

#end of the class definition

# Gets taxonomic info from entrez and adds it to the file formatted in mpa format
def get_mpa_string_From_NCBI(taxID):
	if taxID == 0:
		return "|0:unclassified:U"
	elif taxID == 1:
		return "|1:root:-"
	Entrez.email = getpass.getuser()
	#Creates the data structure from a pull from entrez nucleotide database using accession id with return type of genbank text mode
	handle = Entrez.efetch(db="taxonomy", id=taxID, mode="text", rettype="xml")
	#Parses the returned output into lines
	result= Entrez.read(handle)
	#Will have to change this region to be able to handle more descriptive lists downstream
	recognized_ranks={"cellular organisms":"-", "superkingdom":"d", "kingdom":"k", "phylum":"p", "class":"c", "order":"o", "family":"f", "genus":"g", "species":"s"}
	for entry in result:
		#taxid = entry["Rank"]
		#print(entry)
		special_mpa_string=""
		for r in entry["LineageEx"]:
			if r["Rank"] in recognized_ranks.keys():
				current_rank=recognized_ranks[r["Rank"]]
			else:
				current_rank="-"
			special_mpa_string+=r["TaxId"]+":"+r["ScientificName"]+"|"
	if entry["Rank"] in recognized_ranks.keys():
		current_rank=recognized_ranks[entry["Rank"]]
	else:
		current_rank="-"
	return(special_mpa_string+taxID+":"+entry["ScientificName"])

# Orders the list based on read counts
def order_list(input_kraken, output_list):
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
	top_Node = taxon_Node("root", 0, None, None, 0, "-")
	unclass_Node = taxon_Node("unclassified", 0, None, None, 0, "U")
	for key in sorted(mpa_taxon_counts.keys()):
		key_breakdown=key.split("|")
		for index in range(0,len(key_breakdown)):
			if key_breakdown[index] == "0:unclassified":
				unclass_Node.addCounts(mpa_taxon_counts(key))
			elif index == 0:
				name=key_breakdown[index].split(":")[1]
				current_tax=key_breakdown[index].split(":")[0]
				if top_Node.find_taxID(current_tax) is None:
					temp_Node = taxon_Node(name, mpa_taxon_counts(key), None, None, current_tax, "U")

#def organize_mpas(input_kraken, output_mpa):
def make_node_tree():
	# kraken=open(input_kraken,'r')
	# line=kraken.readline().strip()
	# mpa_dict={}
	# mpa_counts={}
	# counter=0
	# while line != '':
	# 	line_sections = line.split("	")
	# 	contig_id = line_sections[1]
	# 	contig_taxID = line_sections[2]
	# 	if contig_taxID not in mpa_dict.keys():
	# 		print("Adding", contig_taxID)
	# 		mpa_dict[contig_taxID]=get_mpa_string_From_NCBI(contig_taxID)
	# 		mpa_counts[contig_taxID]=1
	# 	else:
	# 		if contig_taxID not in mpa_counts.keys():
	# 			#print("Doesnt exist???", contig_taxID)
	# 			mpa_counts[contig_taxID]=1
	# 		else:
	# 			#print("Incrementing:", contig_taxID)
	# 			mpa_counts[contig_taxID]+=1
	# 	line=kraken.readline().strip()
	# kraken.close()
	headNode = taxon_Node("unclassified", 154, None, None, 0, "u")
	dNode = taxon_Node("Bacteria", 10000, None, None, 2, "d")
	pNode = taxon_Node("Proteobacteria", 9988, None, None, 1224, "p")
	cNode = taxon_Node("Gammaproteobacteria", 9988, None, None, 1236, "c")
	oNode = taxon_Node("Enterobacteriales", 9990, None, None, 91347, "o")
	fNode = taxon_Node("Enterobacteriaceae", 9990, None, None, 543, "f")
	gNode = taxon_Node("Escherichia", 1052, None, None, 516, "g")
	gNode2 = taxon_Node("Klebsiella", 8845, None, None, 570, "g")
	sNode = taxon_Node("coli", 1050, None, None, 562, "s")

	sNode.setParent(gNode)
	gNode.addChild(sNode)
	gNode.setParent(fNode)
	fNode.addChild(gNode)
	gNode2.setParent(fNode)
	fNode.addChild(gNode2)
	fNode.setParent(oNode)
	oNode.addChild(fNode)
	oNode.setParent(cNode)
	cNode.addChild(oNode)
	cNode.setParent(pNode)
	pNode.addChild(cNode)
	pNode.setParent(dNode)
	dNode.addChild(pNode)
	dNode.setParent(headNode)
	headNode.addChild(dNode)
	headNode.print_All()
	print("\n\n\n")
	headNode.find_name("coli").addCounts(1450)
	headNode.print_All()
	total_reads=headNode.getCount()+headNode.getChildren()[0].getCount()
	print("Total reads:", total_reads)


	# mpa_taxon_counts={}
	# #print("mpa_dict length:", len(mpa_dict))
	# for key in mpa_dict.keys():
	# 	#print(key, mpa_dict[key])
	# 	taxons=mpa_dict[key].split("|")
	# 	if taxons is not None:
	# 		#print(taxons)
	# 		for i in range(0, len(taxons)):
	# 			#print(taxons[0:i+1])
	# 			if taxons[i] != "" and taxons[i][0:1] != "-":
	# 				#print(taxons[i])
	# 				if "|".join(taxons[0:i+1]) in mpa_taxon_counts:
	# 					#print("Incrementing", "|".join(taxons[0:i+1]), "from", mpa_taxon_counts["|".join(taxons[0:i+1])], "to",  mpa_taxon_counts["|".join(taxons[0:i+1])]+mpa_counts[key])
	# 					mpa_taxon_counts["|".join(taxons[0:i+1])]+=mpa_counts[key]
	# 				else:
	# 					#print("Creating", "|".join(taxons[0:i+1]), "at 1")
	# 					mpa_taxon_counts["|".join(taxons[0:i+1])]=mpa_counts[key]
	# 	else:
	# 		print("Taxons is none")
	#
	# print("mpa_taxon_counts length:", len(mpa_taxon_counts))
	# for key in sorted(mpa_taxon_counts.keys()):
	# 	print(key, mpa_taxon_counts[key])






	#print("mpa_counts length:", len(mpa_counts))
	#for key in mpa_counts.keys():
	#	print(key, mpa_counts[key])


#get_mpa_string_From_NCBI(470, blank_dick)

args = parseArgs()
order_list(args.input, args.output)
#make_node_tree()
