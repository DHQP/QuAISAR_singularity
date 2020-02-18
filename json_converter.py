#!/usr/bin/env python3

#
# Conversion script from plasmidfinder json output to standard text file
# Usage python ./json_plasmidfinder_converter.py -i plasmidfinder_input_file -o output_file
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

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', help='input file', required=True, dest='input_file')
parser.add_argument('-o', '--output', help='output file', required=True, dest='output_file')
parameters=parser.parse_args()

def convert_json_to_text(infile, outfile):
    with open(infile, 'r') as f:
        hits_dict = json.load(f)

    for hit in hits_dict:
        print(hit['results'])

args=parseArgs()
convert_json_to_text(args.input, args.output)
