#!/bin/sh -l

#$ -o clean_list.out
#$ -e clean_list.err
#$ -N clean_list
#$ -cwd
#$ -q short.q

#Import the config file with shortcuts and settings
if [[ ! -f "./config.sh" ]]; then
	cp ./config_template.sh ./config.sh
fi
. ./config.sh

#
# Description: Script to clean any list file of extra newlines and space
#
# Usage ./clean_list.sh path_to_list_file
#
# Output location: same folder as path_to_list
#
# Modules required: None
#
# v1.0 (10/3/2019)
#
# Created by Nick Vlachos (nvx4@cdc.gov)
#

# Checks for proper argumentation
if [[ $# -eq 0 ]]; then
	echo "No argument supplied to $0, exiting"
	exit 1
elif [[ -z "${1}" ]] || [[ ! -f ${1} ]]; then
	echo "Empty sample name supplied to clean_list.sh or list file does not exist, exiting"
	exit 1
# Gives the user a brief usage and help section if requested with the -h option argument
elif [[ "${1}" = "-h" ]]; then
	echo "Usage is ./clean_list.sh  path_to_list"
	exit 0
fi

# Makes a backup of original list file
cp -f ${1} ${1}.original

# Puts list through dos2unix to convert any windows line returns to unix
dos2unix ${1}.original

tr -d ' \t\r\f' < ${1}.original > ${1}
ex -s +'v/\S/d' -cwq ${1}

# Shows all old lines from list
while IFS= read -r line; do
	echo "O:${line}:"
done < "${1}.original"

# Shows all new lines in list
while IFS= read -r line; do
	echo "N:${line}:"
done < "${1}"

#Script exited gracefully (unless something else inside failed)
exit 0
