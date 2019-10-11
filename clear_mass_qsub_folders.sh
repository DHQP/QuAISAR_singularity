#!/bin/sh -l

#$ -o clear_mqs.out
#$ -e clear_mqs.err
#$ -N clear_mqs
#$ -cwd
#$ -q short.q

#Import the config file with shortcuts and settings
if [[ ! -f "./config.sh" ]]; then
	cp ./config_template.sh ./config.sh
fi
. ./config.sh

#
# Description: Script to clean up script and mass qsub folders of proscripts made while mass submitting many parallele jobs
#
# Usage ./clear_mass_qsub_folders 1|2|3 (1-qsub folders, 2-sharescript outs/errs, 3-Both) folder_containing_script_files_to_be_deleted
#
# Output location: No output created
#
# Modules required: None
#
# v1.0.1 (10/9/2019)
#
# Created by Nick Vlachos (nvx4@cdc.gov)
#

# Number regex to test max concurrent submission parametr
number='^[1-3]+$'

# Checks for proper argumentation
if [[ $# -eq 0 ]]; then
	echo "No argument supplied to $0, exiting"
	exit 1
elif [[ "${1}" = "-h" ]]; then
	echo "Usage is ./clear_mass_qsub_folders.sh  1|2|3 (1-qsub folders, 2-sharescript outs/errs, 3-Both) folder_containing_script_files_to_be_deleted"
	echo "Output by default is downloaded to ${processed}/run_ID and extracted to ${processed}/run_ID/sample_name/FASTQs"
	exit 0
elif [[ -z "${2}" ]]; then
	echo "Empty folder supplied to $0, exiting"
	exit 1
elif ! [[ ${1} =~ $number ]]; then
	echo "Arg1 is not a number or is empty. Please input 1,2, or 3 (1-qsub folders, 2-sharescript outs/errs, 3-Both)...exiting"
	exit 2
fi

# Clears out script folder of all .sh, .err, and .out files and the complete folders within each qsub type folder
if [[ ${1} -eq 1 ]] || [[ ${1} -eq 3 ]]; then
	for folder in ${2}/*
	do
		if [[ -d ${folder} ]]; then
			for folder1 in ${folder}
			do
				echo "Deleting every .sh script in ${folder1}"
				rm ${folder1}/*.sh
				rm ${folder1}/*.out
				rm ${folder1}/*.err

				echo "Deleting every .txt in ${folder1}/complete"
				rm ${folder1}/complete/*.txt
			done
		fi
	done
elif [[ ${1} = "-h" ]]; then
	echo "Will clean out script directory and mass qsub folders depending on parameter. 1 for qsub folder, 2, script folder, or 3 for both"
fi

# Deletes all straggling .err and .out files left in the home shareScript directory
if [[ ${1} -eq 2 ]] || [[ ${1} -eq 3 ]]; then
	rm ${shareScript}/blast16sID_*.out
	rm ${shareScript}/blast16sID_*.err
	rm ${shareScript}/blast16s_*.out
	rm ${shareScript}/blast16s_*.err
	rm ${shareScript}/ani_*.out
	rm ${shareScript}/ani_*.err
	rm ${shareScript}/ANI_*.out
	rm ${shareScript}/ANI_*.err
	rm ${shareScript}/BTQC_*.out
	rm ${shareScript}/BTQC_*.err
	rm ${shareScript}/BUSCO_*.out
	rm ${shareScript}/BUSCO_*.err
	rm ${shareScript}/getFASTQR1_*.out
	rm ${shareScript}/getFASTQR1_*.err
	rm ${shareScript}/getFASTQR2_*.out
	rm ${shareScript}/getFASTQR2_*.err
	rm ${shareScript}/csstn_*.out
	rm ${shareScript}/csstn_*.err
	rm ${shareScript}/csstp_*.out
	rm ${shareScript}/csstp_*.err
	rm ${shareScript}/kraka_*.out
	rm ${shareScript}/kraka_*.err
	rm ${shareScript}/krakr_*.out
	rm ${shareScript}/krakr_*.err
	rm ${shareScript}/gott_*.out
	rm ${shareScript}/gott_*.err
	rm ${shareScript}/mlst_*.out
	rm ${shareScript}/mlst_*.err
	rm ${shareScript}/pFinf_*.out
	rm ${shareScript}/pFinf_*.err
	rm ${shareScript}/pFinp_*.out
	rm ${shareScript}/pFinp_*.err
	rm ${shareScript}/SPAdn_*.out
	rm ${shareScript}/SPAdn_*.err
	rm ${shareScript}/SPAdp_*.out
	rm ${shareScript}/SPAdp_*.err
	rm ${shareScript}/plasFlow_*.out
	rm ${shareScript}/plasFlow_*.err
	rm ${shareScript}/pFlow_*.out
	rm ${shareScript}/pFlow_*.err
	rm ${shareScript}/pFinf_*.out
	rm ${shareScript}/pFinf_*.err
	rm ${shareScript}/PROKK_*.out
	rm ${shareScript}/PROKK_*.err
	rm ${shareScript}/PROKK_*.e*
	rm ${shareScript}/srst2AR_*.out
	rm ${shareScript}/srst2AR_*.err
	rm ${shareScript}/srst2MLST_*.out
	rm ${shareScript}/srst2MLST_*.err
	rm ${shareScript}/srst22MLST_*.out
	rm ${shareScript}/srst22MLST_*.err
	rm ${shareScript}/QUAST_*.out
	rm ${shareScript}/QUAST_*.err
	rm ${shareScript}/QC_*.out
	rm ${shareScript}/QC_*.err
	rm ${shareScript}/MLST_*.out
	rm ${shareScript}/MLST_*.err
	rm ${shareScript}/taxID_*.out
	rm ${shareScript}/taxID_*.err
	rm ${shareScript}/validate_*.out
	rm ${shareScript}/validate_*.err
	rm ${shareScript}/sum_*.out
	rm ${shareScript}/sum_*.err
	rm ${shareScript}/pFn_*.out
	rm ${shareScript}/pFn_*.err
	rm ${shareScript}/pFp_*.out
	rm ${shareScript}/pFp_*.err
	rm ${shareScript}/aniB_*.out
	rm ${shareScript}/aniB_*.err
	rm ${shareScript}/aniM_*.out
	rm ${shareScript}/aniM_*.err
	rm ${shareScript}/node_*.out
	rm ${shareScript}/node_*.err
	rm ${shareScript}/core.*
	rm ${shareScript}/quaisar_*.out
	rm ${shareScript}/quaisar_*.err
	rm ${shareScript}/SNVPhyl_*.out
	rm ${shareScript}/SNVPhyl_*.err
fi
