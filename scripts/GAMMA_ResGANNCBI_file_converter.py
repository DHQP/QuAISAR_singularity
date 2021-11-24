import sys
import glob
import argparse

def ResGANNCBI_Line_Converter(gamma_line):
    Out = ''
    List1 = gamma_line.split('\t')
    Name = List1[0].split('__')
    Out = Name[-1] + '\t' + Name[4] + '\t' + Name[1] + '\t' + Name[2] + '\t'
    Line = '\t'.join(List1[1:])
    Out = Out + Line
    return Out

def ResGANNCBI_File_Converter(gamma):
    f = open(gamma, 'r')
    List1 = []
    for lines in f:
        List1.append(lines)
    f.close()
    Out = open(gamma, 'w')
    Out.write('DB\tResistance\tGene_Family\t' + List1[0])
    for line in List1[1:]:
        New = ResGANNCBI_Line_Converter(line)
        Out.write(New)
    Out.close()

ResGANNCBI_File_Converter(sys.argv[1])
