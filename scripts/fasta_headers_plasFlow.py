#!/usr/bin/env python3

#
# Description: Changes headers in plasFlow assembly fasta from contig# length=length# depth=depthx to Name_contig#_length_length#_depth_depthx
#
# Usage: ./fasta_headers_plasFlow.py -i input.fasta -o output.fasta
#
# Output location: parameter
#
# Modules required: Biopython must be available in python instance
#
# V1.0
#
# Created by Erisa Sula (nvd4@cdc.gov) and Nick Vlachos (nvx4@cdc.gov)
#

from Bio import SeqIO
import sys
import os
import argparse

#Create an arg parser...someday
def parseArgs(args=None):
	parser = argparse.ArgumentParser(description='Script to rename contigs in re-assembled plasFlow filtered assemblies')
	parser.add_argument('-i', '--input', required=True, help='input fasta filename')
	parser.add_argument('-o', '--output', required=True, help='output filename')
	return parser.parse_args()

args = parseArgs()
sequences = []
for record in SeqIO.parse(args.input,"fasta"):
    #print(record.id)
    name=os.path.basename(args.input).split("_")[::-1]
    name=name[2:]
    name='_'.join(name[::-1])
    #print(name)
    #record.id = record.id.split("_cov")[0].replace("NODE",name)
    print(record)
    print(name)
    print(record.description.split(" ")[0])
    contig = record.description.split(" ")[0]
    record.description.split(" ")[1]
    length = record.description.split(" ")[1].split("=")[1]
    depth = record.description.split(" ")[2].split("=")[1]
    record.id = name+"_"+contig+"_length_"+length+"_depth_"+depth

    #print(record.id)
    record.description = ""
#    print(record.description)
#    print(record)
    sequences.append(record)

SeqIO.write(sequences, args.output, "fasta")
