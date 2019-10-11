#!/usr/bin/env python3

#
# Description: Script to choose the centroid sample, via mash distances, within a list of samples
#
# Usage: python3 ./Mash_centroid.py -i input_list_file -o output_list_filename
#
# Output location: parameter
#
# Modules required: None
#
# v1.0.1 (10/7/2019)
#
# Created by Rich Stanton (njr5@cdc.gov)
#

import sys
import glob
import numpy
import operator
from operator import itemgetter
import subprocess
import argparse

# Parse all argument from command line
def parseArgs(args=None):
	parser = argparse.ArgumentParser(description='Script to choose the centroid sample, via mash distances, within a list of samples')
	parser.add_argument('-i', '--input', required=True, help='input list')
	parser.add_argument('-o', '--output', required=True, help='output centroided lst')
	return parser.parse_args()

def Mash_List(Mash_Index):
    """Takes in an index of Mash files and makes a list of the info from each"""
    List1 = glob.glob(Mash_Index)
    Combined = []
    for files in List1:
        f = open(files, 'r')
        String1 = f.readline()
        Combined.append(String1)
        f.close()
    Combined.sort()
    return Combined

def Average_Mash(input_mash_list):
    """Takes in a Mash_List (made using Mash_List) and then makes a list with the average value for each entry"""
    Average_List = []
    List1 = input_mash_list[0].split('\t')
    Entry = List1[0]
    values = []
    for entries in input_mash_list:
        List1 = entries.split('\t')
        if List1[0] != Entry:
            average = numpy.mean(values)
            Total = []
            Total.append(Entry)
            Total.append(average)
            Average_List.append(Total)
            Entry = List1[0]
            values = []
            values.append(float(List1[2]))
        else:
            values.append(float(List1[2]))
    average = numpy.mean(values)
    Total = []
    Total.append(Entry)
    Total.append(average)
    Average_List.append(Total)
    Average_List.sort(key=operator.itemgetter(1))
    return Average_List

def Mash_List_Maker(input_assembly_list):
    """Makes a list of tha all x all mash outputs"""
    Output = []
    for files in input_assembly_list:
        for files2 in input_assembly_list:
            String1 = subprocess.check_output('mash dist ' + files + ' ' + files2, shell=True)
            Output.append(String1)
    Output.sort()
    return Output

def Mash_Centroid(input_assembly_list):
    """Returns the name of the fasta with the lowest average mash index"""
    List1 = Mash_List_Maker(input_assembly_list)
    Averages = Average_Mash(List1)
    Best = Averages[0][0]
    return Best

def Fasta_List(input_file):
    """Takes in an input file for scicomp isolates and makes a list of their assembly_locations"""
    Out_List = []
    f = open(input_file, 'r')
    String1 = f.readline()
    while String1 != '':
        if String1[-1] == '\n':
            String1 = String1[0:-1]
        List1 = String1.split('/')
        Out_String = '/scicomp/groups/OID/NCEZID/DHQP/CEMB/MiSeqAnalysisFiles/' + String1 + '/Assembly/' + List1[1] + '_scaffolds_trimmed.fasta'
        Out_List.append(Out_String)
        String1 = f.readline()
    f.close()
    return Out_List

def List_Reorder(input_file, best_isolate, output_file):
    """Takes in an input_list file and a best fasta and makes an output file with the best isolate on top"""
    Output = open(output_file, 'w')
    f = open(input_file, 'r')
    Output.write(best_isolate + '\n')
    String1 = f.readline()
    while String1 != '':
        if (best_isolate in String1) == True:
            String1 = f.readline()
        else:
            Output.write(String1)
            String1 = f.readline()
    f.close()
    Output.close()

def Scicomp_Mash_Centroid(input_list, output_list):
    """Takes in an input list and returns an output list with the centroid isolate at the top"""
    Fastas = Fasta_List(input_list)
    Centroid = Mash_Centroid(Fastas)
    Best_List = Centroid.split('/')
    Best_Centroid = Best_List[8] + '/' + Best_List[9]
    List_Reorder(input_list, Best_Centroid, output_list)

args = parseArgs()
Scicomp_Mash_Centroid(args.input, args.output)
