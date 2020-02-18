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
    with open(infile, 'r') as file:
        data = file.read().replace('{"plasmidfinder": {"results": ', '')[:-2]
    results_json=json.loads(data)



#     with open(infile, 'r') as f:
#         top_layer = json.load(f)
#
#     for keys,values in top_layer.items():
#         print(keys)
#         print(values)

    for keys,values in results_json.items():
        print(keys)
        print(values)


args = parseArgs()
convert_json_to_text(args.input, args.output)
