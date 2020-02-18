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
import argparse

def parseArgs(args=None):
    parser = argparse.ArgumentParser(description='Script to parse through plasmidfinder results')
    parser.add_argument('-i', '--input', help='input file', required=True)
    parser.add_argument('-o', '--output', help='output file', required=True)
    return parser.parse_args()

def convert_json_to_text(infile, outfile):
    with open(infile, 'r') as f:
        results = json.load(f)

    print(results['plasmidfinder']['results']['Gram Positive'])

    for key,value in results['plasmidfinder']['results']['Enterobacteriaceae']['enterobacteriaceae'].items():
        print(key)
        print(value)

args = parseArgs()
convert_json_to_text(args.input, args.output)
