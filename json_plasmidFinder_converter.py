#!/usr/bin/env python3

#
# Conversion script from plasmidFinder json output to standard text file
# Usage python ./json_plasmidFinder_converter.py -i plasmidFinder_input_file -o output_file
#
# Output location: Parameter
#
# Modules required: None
#
# v1.0 (02/18/2020)
#
# Created by Nick Vlachos (nvx4@cdc.gov)
#

import json
import argparse

def parseArgs(args=None):
    parser = argparse.ArgumentParser(description='Script to parse through plasmidFinder results')
    parser.add_argument('-i', '--input', help='input file', required=True)
    parser.add_argument('-o', '--output', help='output file', required=True)
    return parser.parse_args()

def convert_json_to_text(infile, outfile):
    gramps=[]
    entero=[]

    with open(infile, 'r') as f:
        results = json.load(f)

    # Parse all Gram positive hits
    for key,value in results['plasmidfinder']['results']['Gram Positive'].items():
        #print(key)
        #print(value)
        if value != "No hit found":
            plasmid = value['plasmid']
            percent_identity = str(value['identity'])
            HSP_length = str(value['HSP_length'])
            template_length = str(value['template_length'])
            contig = value['contig_name']
            contig_position = value['positions_in_contig']
            accession_number = value['note']
            coverage = str(value['coverage'])
            gramps.append(plasmid+"\t"+percent_identity+"\t"+HSP_length+'/'+template_length+"\t"+contig+"\t"+contig_position+"\t"+accession_number+"\n")


    # Parse all Enterobacteriaceae hits
    print(results['plasmidfinder']['results']['Enterobacteriaceae']['enterobacteriaceae'])
    for key,value in results['plasmidfinder']['results']['Enterobacteriaceae']['enterobacteriaceae'].items():
        #print(key)
        #print(value)
        if value != "No hit found":
            plasmid = value['plasmid']
            percent_identity = str(value['identity'])
            HSP_length = str(value['HSP_length'])
            template_length = str(value['template_length'])
            contig = value['contig_name']
            contig_position = value['positions_in_contig']
            accession_number = value['note']
            coverage = str(value['coverage'])
            entero.append(plasmid+"\t"+percent_identity+"\t"+HSP_length+'/'+template_length+"\t"+contig+"\t"+contig_position+"\t"+accession_number+"\n")
    # Add nothing found if no hits
    if len(gramps) == 0:
        gramps.append("No plasmid replicons found.\n")
    if len(entero) == 0:
        entero.append("No plasmid replicons found.\n")

    # print("::S::")
    # for line in gramps:
    #     print(line)
    # for line in entero:
    #     print(line)
    # print("::F::")

    # Write findings to file
    f=open(outfile, "w")
    f.write("Plasmid\tidentity\tHSP/Query\tcontig\tPosition_in_contig\tAccesion\n")
    f.write("Enterococcus,Streptococcus,Staphylococcus\n")
    for line in gramps:
        f.write(line)
    f.write("Enterobacteriaceae\n")
    for line in entero:
        f.write(line)
    f.close()

args = parseArgs()
convert_json_to_text(args.input, args.output)
