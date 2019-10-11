#!/bin/python3

#
# Description: A script that takes a sample and compares it to a busco database to discover number of similar genes (% conserved proteins) from prokka output
#
# Usage python3 ./do_busco.py -n sample_name -d DB(to search against) -p run_ID
#
# Output location: default_config.sh_output_location/run_ID/sample_name/BUSCO/
#
# Modules required: None
#
# v1.0 (10/9/2019)
#
# Created by Nick Vlachos (nvx4@cdc.gov)
#

import os
import sys
import argparse
import shutil
import subprocess


#Create an arg parser...someday
def parseArgs(args=None):
	parser = argparse.ArgumentParser(description='Tool to call BUSCO on a sample')
	parser.add_argument('-n', '--name', required=True, help='Name of isolate', dest="name")
	parser.add_argument('-d', '--database', required=True, help='database to run against', dest="db")
	parser.add_argument('-p', '--project', required=True, help='project ID', dest="project")
	return parser.parse_args()

args = parseArgs()

# Sets output folder to parent isolate folder
OUTDATADIR="/scicomp/groups/OID/NCEZID/DHQP/CEMB/MiSeqAnalysisFiles/"+args.project+"/"+args.name+"/"
# Sets buscoDB to the folder matching the 2nd argument (dbs vary on taxonomic level)
buscoDB="/scicomp/groups/OID/NCEZID/DHQP/CEMB/BUSCO/"+args.db.lower()+"/"

### BUSCO % Conserved Protein Identity ###
print("Running BUSCO for Conserved Protein Identity")

#Deletes old busco folder if existent
if os.path.isdir(OUTDATADIR+"BUSCO"):
	shutil.rmtree(OUTDATADIR+"BUSCO")


# Creates busco output folder if not already existent
print ("Creating "+OUTDATADIR+"BUSCO")
os.mkdir(OUTDATADIR+"BUSCO")

# Get current directory, as it needs to be changed to control busco output
owd=os.getcwd()
os.chdir(OUTDATADIR)

# Run busco on the prokka output using the database provided by command line ($2)
print(sys.version)

#prokka = OUTDATADIR+"prokka/"+args.name+"_PROKKA.faa"
subprocess.call(["run_BUSCO.py", "-i", OUTDATADIR+"prokka/"+args.name+"_PROKKA.faa", "-o", args.name, "-l", args.db, "-m", "prot"], shell=True)

# Moves output files to proper location and removes temp files
source = OUTDATADIR+"run_"+args.name+"/"
dest = OUTDATADIR+"BUSCO/"
files = os.listdir(source)
for f in files:
        shutil.move(source, dest)

#Deletes temp folder_containing_script_files_to_be_deleted
shutil.rmtree(OUTDATADIR+"run_"+args.name)
shutil.rmtree(OUTDATADIR+"tmp")

# returns current directory to original location
os.chdir(owd)
