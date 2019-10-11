#!/usr/bin/env python3

#
# Description: Script to summarize the weighted kraken file output
#
# Usage: python3 ./Kraken_Assembly_Summary_Exe.py -k input_kraken_file -l label_file -t list_file -o output_filename
#
# Output location: parameter
#
# Modules required: Biopython must be available in python instance
#
# v1.0 (10/3/2019)
#
# Created by Rich Stanton (njr5@cdc.gov)
#

import sys
import glob
from Bio import SeqIO
import argparse
from decimal import Decimal, ROUND_HALF_UP

def parseArgs(args=None):
	parser = argparse.ArgumentParser(description='Script to summarize the kraken output files for weighted versions')
	parser.add_argument('-k', '--kraken', required=True, help='input weighted kraken filename')
	parser.add_argument('-l', '--label', required=True, help='input label filename')
	parser.add_argument('-t', '--list', required=True, help='input list filename')
	parser.add_argument('-o', '--output', required=True, help='output list filename')
	return parser.parse_args()

def In_List(item, list1):
    """determines if item is in list1"""
    for x in list1:
        if x == item:
            return True
    else: return False

def In_Line(value, line_string):
    if len(str(value)) > len(line_string):
        return False
    else:
        for characters in range(0, (len(line_string)-len(value) + 1)):
            if value == line_string[characters:(characters + len(value))]:
                return True
        return False

def String_Converter(Input_String):
    Counter = 0
    Character = Input_String[0]
    Current_String = ''
    Out_List = []
    while Counter < len(Input_String):
        if Character == '\n':
            Out_List.append(Current_String)
            return Out_List
        elif Character == '\t':
            Out_List.append(Current_String)
            Current_String = ''
            Counter = Counter + 1
            Character = Input_String[Counter]
        else:
            Current_String = Current_String + Character
            Counter = Counter + 1
            Character = Input_String[Counter]

def Translated_Line_Converter(Input_String):
    Counter = 0
    Character = Input_String[0]
    Current_String = ''
    Out_List = []
    while Counter < len(Input_String):
        if Character == '\n':
            Out_List.append(Current_String)
            return Out_List
        elif Character == '\t' or Character == ';':
            Out_List.append(Current_String)
            Current_String = ''
            Counter = Counter + 1
            Character = Input_String[Counter]
        else:
            Current_String = Current_String + Character
            Counter = Counter + 1
            Character = Input_String[Counter]

def Minimum_Contig(input_kraken):
    f = open(input_kraken, 'r')
    String1 = f.readline()
    List1 = String_Converter(String1)
    Minimum = int(List1[3])
    while String1 != '':
        List1 = String_Converter(String1)
        if Minimum > int(List1[3]):
            Minimum = int(List1[3])
            String1 = f.readline()
        else:
            String1 = f.readline()
    f.close()
    return Minimum

def Kraken_Assembly_Converter(input_kraken, output_kraken):
    """Reads in a kraken input file and an input fasta and makes a bp weighted kraken file"""
    Min_Length = Minimum_Contig(input_kraken)
    f = open(input_kraken, 'r')
    g = open(output_kraken, 'w')
    String1 = f.readline()
    while String1 != '':
        List1 = String_Converter(String1)
        Count = int(round(float(List1[3]) // Min_Length))
        Out_Line = ''
        for items in List1[0:4]:
            Out_Line = Out_Line + items + '\t'
        Out_Line = Out_Line[0:-1] + '\n'
        for numbers in range(Count):
            g.write(Out_Line)
        String1 = f.readline()
    f.close()
    g.close()

def Summary_Taxonomy(input_line):
    Output = ''
    for characters in input_line:
        if characters == '\n':
            return Output
        elif characters == '\t':
            Output = ''
        elif characters == ' ' and Output == '':
            Output = ''
        else:
            Output = Output + characters

def Kraken_Assembly_Lengths(input_kraken):
    Output_List = []
    f = open(input_kraken, 'r')
    String1 = f.readline()
    while String1 != '':
        List1 = String_Converter(String1)
        Entry = []
        Entry.append(List1[1])
        Entry.append(List1[3])
        Output_List.append(Entry)
        String1 = f.readline()
    f.close()
    return Output_List

def Kraken_Assembly_Unclassified(input_kraken):
    Length = 0
    f = open(input_kraken)
    String1 = f.readline()
    while String1 != '':
        List1 = String_Converter(String1)
        if List1[2] == '0':
            Length = Length + int(List1[3])
        String1 = f.readline()
    f.close()
    return Length

def Total_Length(input_kraken):
    f = open(input_kraken)
    Length = 0
    String1 = f.readline()
    while String1 != '':
        List1 = String_Converter(String1)
        Length = Length + int(List1[3])
        String1 = f.readline()
    f.close()
    return Length



def Taxonomy_Totals(untranslated_kraken, translated_kraken):
    Contig_Lengths = Kraken_Assembly_Lengths(untranslated_kraken)
    f = open(translated_kraken, 'r')
    Entry_List = []
    String1 = f.readline()
    while String1 != '':
        List1 = Translated_Line_Converter(String1)
        for contigs in Contig_Lengths:
            if contigs[0] == List1[0]:
                for classes in range(1, (len(List1) - 1)):
                    Class = []
                    Class.append(List1[classes])
                    Class.append(contigs[1])
                    Class.append(0)
                    Entry_List.append(Class)
                Class = []
                Class.append(List1[len(List1) - 1])
                Class.append(contigs[1])
                Class.append(contigs[1])
                Entry_List.append(Class)
            else:
                continue
        String1 = f.readline()
    f.close()
    Unclassified_Length = Kraken_Assembly_Unclassified(untranslated_kraken)
    Unclassified = ['unclassified']
    Unclassified.append(str(Unclassified_Length))
    Unclassified.append(str(Unclassified_Length))
    Entry_List.append(Unclassified)
    return Entry_List

def List_Combiner(input_list):
    Output_List = []
    Output_List.append(input_list[0])
    for items in range(1, len(input_list)):
        Match = 0
        for entries in Output_List:
            if input_list[items][0] == entries[0]:
                entries[1] = str(int(entries[1]) + int(input_list[items][1]))
                entries[2] = str(int(entries[2]) + int(input_list[items][2]))
                Match = 1
        if Match == 0:
            Output_List.append(input_list[items])
    return Output_List

def Two_Decimal_Percent(number, total):
    our_value = 100 * Decimal(float(number) // total)
    output = Decimal(our_value.quantize(Decimal('0.01'), rounding=ROUND_HALF_UP))
    Value = str(output)
    return Value

def Assembly_Summary(input_summary, input_list, output_summary):
    Total_Length = 0
    for items in input_list:
        Total_Length = Total_Length + int(items[2])
    f= open(input_summary, 'r')
    g = open(output_summary, 'w')
    String1 = f.readline()
    while String1 != '':
        List1 = String_Converter(String1)
        for entries in input_list:
            if entries[0] == Summary_Taxonomy(String1):
                Out_Line = ''
                Percent = Two_Decimal_Percent(int(entries[1]), Total_Length)
                if len(Percent) == 4:
                    Out_Line = '  ' + Percent + '\t' + entries[1] + '\t' + str(entries[2]) + '\t' + List1[3] + '\t' + List1[4] + '\t' + List1[5] + '\n'
                elif len(Percent) == 5:
                    Out_Line = ' ' + Percent + '\t' + entries[1] + '\t' + str(entries[2]) + '\t'+ List1[3] + '\t' + List1[4] + '\t' + List1[5] + '\n'
                elif len(Percent) == 6:
                    Out_Line = Percent + '\t' + entries[1] + '\t' + str(entries[2]) + '\t'+ List1[3] + '\t' + List1[4] + '\t' + List1[5] + '\n'
                g.write(Out_Line)
                break
        String1 = f.readline()
    f.close()
    g.close()

args = parseArgs()
List2 = Taxonomy_Totals(args.kraken, args.label)
List3 = List_Combiner(List2)
Assembly_Summary(args.list, List3, args.output)


##Kraken_Assembly_Converter(sys.argv[1], sys.argv[2], sys.argv[1][0:-3] + '_BP.in')
