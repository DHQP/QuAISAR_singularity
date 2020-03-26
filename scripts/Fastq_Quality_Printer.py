#!/usr/bin/env python3

#
# Description: Changes headers in SPAdes assembly fasta from contig# length=length# depth=depthx to Name_contig#_length_length#_depth_depthx
#
# Output location: standard out
#
# Usage: ./Fastq_Quality_Printer.py -1 path_to_R1 -2 path_to_R2
#
# Modules required: None
#
# v1.0 (10/3/2019)
#
# Created by Rich Stanton (njr5@cdc.gov)
#

import sys
import argparse
from decimal import *
getcontext().prec = 4

def parseArgs(args=None):
	parser = argparse.ArgumentParser(description='Script to count quality metrics from paired fastq files')
	parser.add_argument('-1', '--r1', required=True, help='input R1 filename')
	parser.add_argument('-2', '--r2', required=True, help='input R2 filename')
	return parser.parse_args()


def Q20(input_string):
    Q20_Total = 0
    for letters in input_string:
        if (ord(letters) - 33) >=20:
            Q20_Total = Q20_Total + 1
    return Q20_Total

def Q30(input_string):
    Q30_Total = 0
    for letters in input_string:
        if (ord(letters) - 33) >=30:
            Q30_Total = Q30_Total + 1
    return Q30_Total

def Quality_Score(input_fastq):
    f= open(input_fastq, 'r')
    Q20_Total = 0
    Q30_Total = 0
    Total_Bases = 0
    Total_Reads = 0
    String1 = f.readline()
    while String1 != '':
        if String1 == '+\n':
            Total_Reads = Total_Reads + 1
            String1 = f.readline()
            Q20_Total = Q20_Total + Q20(String1[0:-1])
            Q30_Total = Q30_Total + Q30(String1[0:-1])
            Total_Bases = Total_Bases + len(String1[0:-1])
            String1 = f.readline()
        else:
            String1 = f.readline()
    print('Total Reads: ' + str(Total_Reads))
    print('Total Bases: ' + str(Total_Bases))
    print('Q20 Bases: ' + str(Q20_Total))
    print('Q30 Bases: ' + str(Q30_Total))
    print('Q20 %: ' + str(float(Q20_Total) / Total_Bases))
    print('Q30 %: ' + str(float(Q30_Total) / Total_Bases))

def Quality_Score_2_Reads(fastq1, fastq2, output_file):
    Q20_Total_1 = 0
    Q30_Total_1 = 0
    Total_Bases_1 = 0
    Total_Reads_1 = 0
    f = open(fastq1, 'r')
    String1 = f.readline()
    while String1 != '':
        if String1 == '+\n':
            Total_Reads_1 = Total_Reads_1 + 1
            String1 = f.readline()
            Q20_Total_1 = Q20_Total_1 + Q20(String1[0:-1])
            Q30_Total_1 = Q30_Total_1 + Q30(String1[0:-1])
            Total_Bases_1 = Total_Bases_1 + len(String1[0:-1])
            String1 = f.readline()
        else:
            String1 = f.readline()
    f.close()
    Q20_Total_2 = 0
    Q30_Total_2 = 0
    Total_Bases_2 = 0
    Total_Reads_2 = 0
    g = open(fastq2, 'r')
    String1 = g.readline()
    while String1 != '':
        if String1 == '+\n':
            Total_Reads_2 = Total_Reads_2 + 1
            String1 = g.readline()
            Q20_Total_2 = Q20_Total_2 + Q20(String1[0:-1])
            Q30_Total_2 = Q30_Total_2 + Q30(String1[0:-1])
            Total_Bases_2 = Total_Bases_2 + len(String1[0:-1])
            String1 = g.readline()
        else:
            String1 = g.readline()
    g.close()
    Total_Reads = str(Total_Reads_1 + Total_Reads_2)
    Total_Bases = str(Total_Bases_1 + Total_Bases_2)
    Q20_Total = str(Q20_Total_1 + Q20_Total_2)
    Q30_Total = str(Q30_Total_1 + Q30_Total_2)
    Q20_R1 = str(Decimal(Q20_Total_1) / Decimal(Total_Bases_1))
    Q20_R2 = str(Decimal(Q20_Total_2) / Decimal(Total_Bases_2))
    Q30_R1 = str(Decimal(Q30_Total_1) / Decimal(Total_Bases_1))
    Q30_R2 = str(Decimal(Q30_Total_2) / Decimal(Total_Bases_2))
    h = open(output_file, 'w')
    h.write(fastq1 + '\t' + Q20_Total + '\t' + Q30_Total + '\t' + str(Q20_Total_1) + '\t' + str(Q20_Total_2) + '\t' + Q20_R1 + '\t' + Q20_R2 + '\t' + str(Q30_Total_1) + '\t' + str(Q30_Total_2) + '\t' + Q30_R1 + '\t' + Q30_R2+ '\t' + Total_Bases + '\t' + Total_Reads)
    h.close()

def Quality_Score_Printer(fastq1, fastq2):
    Q20_Total_1 = 0
    Q30_Total_1 = 0
    Total_Bases_1 = 0
    Total_Reads_1 = 0
    f = open(fastq1, 'r')
    String1 = f.readline()
    while String1 != '':
        if String1.startswith('+'):
            Total_Reads_1 = Total_Reads_1 + 1
            String1 = f.readline()
            Q20_Total_1 = Q20_Total_1 + Q20(String1[0:-1])
            Q30_Total_1 = Q30_Total_1 + Q30(String1[0:-1])
            Total_Bases_1 = Total_Bases_1 + len(String1[0:-1])
            String1 = f.readline()
        else:
            String1 = f.readline()
    f.close()
    Q20_Total_2 = 0
    Q30_Total_2 = 0
    Total_Bases_2 = 0
    Total_Reads_2 = 0
    g = open(fastq2, 'r')
    String1 = g.readline()
    while String1 != '':
        if String1.startswith('+'):
            Total_Reads_2 = Total_Reads_2 + 1
            String1 = g.readline()
            Q20_Total_2 = Q20_Total_2 + Q20(String1[0:-1])
            Q30_Total_2 = Q30_Total_2 + Q30(String1[0:-1])
            Total_Bases_2 = Total_Bases_2 + len(String1[0:-1])
            String1 = g.readline()
        else:
            String1 = g.readline()
    g.close()
    Total_Reads = str(Total_Reads_1 + Total_Reads_2)
    Total_Bases = str(Total_Bases_1 + Total_Bases_2)
    Q20_Total = str(Q20_Total_1 + Q20_Total_2)
    Q30_Total = str(Q30_Total_1 + Q30_Total_2)
    Q20_R1 = str(Decimal(Q20_Total_1) / Decimal(Total_Bases_1))
    Q20_R2 = str(Decimal(Q20_Total_2) / Decimal(Total_Bases_2))
    Q30_R1 = str(Decimal(Q30_Total_1) / Decimal(Total_Bases_1))
    Q30_R2 = str(Decimal(Q30_Total_2) / Decimal(Total_Bases_2))
    String1 = fastq1 + '\t' + Q20_Total + '\t' + Q30_Total + '\t' + str(Q20_Total_1) + '\t' + str(Q20_Total_2) + '\t' + Q20_R1 + '\t' + Q20_R2 + '\t' + str(Q30_Total_1) + '\t' + str(Q30_Total_2) + '\t' + Q30_R1 + '\t' + Q30_R2+ '\t' + Total_Bases + '\t' + Total_Reads
    print(String1)


args = parseArgs()
Quality_Score_Printer(args.r1, args.r2)
