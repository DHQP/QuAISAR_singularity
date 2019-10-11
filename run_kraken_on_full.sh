#!/bin/sh -l

#$ -o run_kraken.out
#$ -e run_kraken.err
#$ -N run_kraken
#$ -cwd
#$ -q short.q

#Import the config file with shortcuts and settings
if [[ ! -f "./config.sh" ]]; then
	cp ./config_template.sh ./config.sh
fi
. ./config.sh

#
# Description: Runs a fasta (assembly) through the kraken classification tool which identifies the most likely
# 	taxonomic classification for the sample on the full kraken database
#
# Usage: ./run_kraken_on_full.sh sample_name run_ID
#
# Output location: default_config.sh_output_location/run_ID/sample_name/kraken_full/
#
# Modules required: kraken/0.10.5, perl/5.12.3, krona/2.7, Python3/3.5.2
#
# v1.0 (10/3/2019)
#
# Created by Nick Vlachos (nvx4@cdc.gov)
#

ml kraken/0.10.5 perl/5.12.3 krona/2.7 Python3/3.5.2

# Checks for proper argumentation
if [[ $# -eq 0 ]]; then
	echo "No argument supplied to $0, exiting"
	exit 1
elif [[ -z "${1}" ]]; then
	echo "Empty sample name supplied to run_kraken.sh, exiting"
	exit 1
# Gives the user a brief usage and help section if requested with the -h option argument
elif [[ "${1}" = "-h" ]]; then
	echo "Usage is ./run_kraken.sh   sample_name	run_ID"
	echo "Output is saved to in ${processed}/miseq_run_ID/sample_name/kraken/full"
	exit 0
elif [ -z "$2" ]; then
	echo "Empty assembly relativity supplied to run_kraken.sh. Second argument should be 'pre' for paired reads or 'post' for assembly (no quotes). Exiting"
	exit 1
fi

# Sets the parent output folder as the sample name folder in the processed samples folder in MMB_Data
OUTDATADIR="${processed}/${4}/${1}"

# Creates folder for output from kraken
if [ ! -d "$OUTDATADIR/kraken" ]; then
	echo "Creating $OUTDATADIR/kraken"
	mkdir -p "$OUTDATADIR/kraken/postAssembly_full"
elif [ ! -d "$OUTDATADIR/kraken/postAssembly_full" ]; then
	echo "Creating $OUTDATADIR/kraken/postAssembly_full"
	mkdir -p "$OUTDATADIR/kraken/postAssembly_full"
fi

# Prints out version of kraken
kraken --version
# Status view of run
echo "[:] Running kraken.  Output: ${1}.kraken / ${1}.classified"
# Runs chosen kraken db on the assembly
kraken --db "${kraken_full_db}" --preload --threads "${procs}" --output "${OUTDATADIR}/kraken/postAssembly_full/${1}_full.kraken" --classified-out "${OUTDATADIR}/kraken/postAssembly_full/${1}_full.classified" "${OUTDATADIR}/Assembly/${1}_scaffolds_trimmed.fasta"
# Attempting to weigh contigs and produce standard krona and list output using a modified version of Rich's weighting scripts (will also be done on pure contigs later)
python ${shareScript}/Kraken_Assembly_Converter_2_Exe.py -i "${OUTDATADIR}/kraken/postAssembly_full/${1}_full.kraken"
kraken-translate --db "${kraken_full_db}" "${OUTDATADIR}/kraken/postAssembly_full/${1}_full_BP.kraken" > "${OUTDATADIR}/kraken/postAssembly_full/${1}_full_BP.labels"
kraken-mpa-report --db "${kraken_full_db}" "${OUTDATADIR}/kraken/postAssembly_full/${1}_full_BP.kraken" > "${OUTDATADIR}/kraken/postAssembly_full/${1}_full_weighted.mpa"
python3 "${shareScript}/Metaphlan2krona.py" -p "${OUTDATADIR}/kraken/postAssembly_full/${1}_full_weighted.mpa" -k "${OUTDATADIR}/kraken/postAssembly_full/${1}_full_weighted.krona"
kraken-report --db "${kraken_full_db}" "${OUTDATADIR}/kraken/postAssembly_full/${1}_full_BP.kraken" > "${OUTDATADIR}/kraken/postAssembly_full/${1}_full_BP.list"
python ${shareScript}/Kraken_Assembly_Summary_Exe.py -k "${OUTDATADIR}/kraken/postAssembly_full/${1}_full_BP.kraken" -l "${OUTDATADIR}/kraken/postAssembly_full/${1}_full_BP.labels" -t "${OUTDATADIR}/kraken/postAssembly_full/${1}_full_BP.list" -o "${OUTDATADIR}/kraken/postAssembly_full/${1}_full_BP_data.list"
#. "${shareScript}/module_changers/perl_5221_to_5123.sh"
ktImportText "${OUTDATADIR}/kraken/postAssembly_full/${1}_full_weighted.krona" -o "${OUTDATADIR}/kraken/postAssembly_full/${1}_full_weighted_BP_krona.html"
#. "${shareScript}/module_changers/perl_5123_to_5221.sh"
"${shareScript}/best_hit_from_kraken.sh" "${1}" "${2}" "full_BP_data" "${4}"

# Run the metaphlan generator on the kraken output
echo "[:] Generating metaphlan compatible report."
kraken-mpa-report --db "${kraken_full_db}" "${OUTDATADIR}/kraken/postAssembly_full/${1}_full.kraken" > "${OUTDATADIR}/kraken/postAssembly_full/${1}_full.mpa"

# Run the krona generator on the metaphlan output
echo "[:] Generating krona output for ${1}."
python3 "${shareScript}/Metaphlan2krona.py" -p "${OUTDATADIR}/kraken/postAssembly_full/${1}_full.mpa" -k "${OUTDATADIR}/kraken/postAssembly_full/${1}_full.krona"

# Run the krona graph generator from krona output
#. "${shareScript}/module_changers/perl_5221_to_5123.sh"
ktImportText "${OUTDATADIR}/kraken/postAssembly_full/${1}_full.krona" -o "${OUTDATADIR}/kraken/postAssembly_full/${1}_full.html"
#. "${shareScript}/module_changers/perl_5123_to_5221.sh"

# Creates the parsible report from the kraken output
echo "[:] Creating alternate report for taxonomic extraction"
kraken-report --db "${kraken_full_db}" "${OUTDATADIR}/kraken/postAssembly_full/${1}_full.kraken" > "${OUTDATADIR}/kraken/postAssembly_full/${1}_full.list"

# Parses the output for the best taxonomic hit
echo "[:] Extracting best taxonomic matches"

# Runs the extractor for pulling best hit from a kraken run
"${shareScript}/best_hit_from_kraken.sh" "${1}" "${2}" "full" "${4}" "${5}"

ml -kraken/0.10.5 -perl/5.12.3 -krona/2.7 -Python3/3.5.2

#Script exited gracefully (unless something else inside failed)
exit 0
