#! /usr/bin/env bash
#!/bin/sh -l


#
# Description: Script to prepare environment to be able to run quaisar pipeline
#
# Usage: ./quaisar_installation.sh location_for
#
# Output location: No output created
#
# Modules required: None
#
# v1.0 (3/24/2020)
#
# Created by Nick Vlachos (nvx4@cdc.gov)
#

#finds where script is it, so it properly reference directories during install
install_script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd | rev | cut -d'/' -f2- | rev)"

echo "${install_script_dir}"

# Checking for proper number of arguments from command line
if [[ $# -lt 1  || $# -gt 3 ]]; then
  echo "Usage: ./installation.sh script_installation_location database_installation_location Working_directory_of_output_from_runs"
	echo "You have used $# args"
  exit 3
fi

# Checks for proper argumentation
if [[ $# -eq 0 ]]; then
	echo "No argument supplied to $0, exiting"
	exit 1
elif [ -z "$1" ]; then
	echo "Empty script location supplied to supplied to $0, exiting"
	exit 1
elif [ -z "$2" ]; then
	echo "Empty database installation location supplied to $0, exiting"
	exit 1
elif [ -z "$3" ]; then
  echo "Empty working directory for output supplied to $0, exiting"
  exit 1
# command line version of usage for script
elif [[ "$1" = "-h" ]]; then
	echo "Usage: ./installation.sh script_installation_location database_installation_location Working_directory_of_output_from_runs"
	exit 0
fi



# need to add to bottom of yaml with proper home location
installation_location=${1}
databases=${2}
working_directory=${3}

echo "Installing quaisar scripts to ${installation_location}"
if [[ ! -d ${installation_location} ]]; then
  echo "Creating ${installation}"
  mkdir -p ${installation_location}
fi

echo "Checking for working directory of runs ${working_directory}"
if [[ ! -d ${working_directory} ]]; then
  echo "Creating ${working_directory}"
  mkdir -p ${working_directory}
fi

# Create primary config file
head -n25 ${install_script_dir}/installation/config_template.sh > ${installation_location}/new_config.sh
echo "output_dir=${working_directory}" >> ${installation_location}/new_config.sh
head -n27 ${install_script_dir}/installation/config_template.sh | tail -n1 >> ${installation_location}/new_config.sh
echo "src=${installation_location}" >> ${installation_location}/config.sh
head -n29 ${install_script_dir}/installation/config_template.sh | tail -n1 >> ${installation_location}/new_config.sh
echo "local_DBs=${databases}" >> ${installation_location}/config.sh
head -n31 ${install_script_dir}/installation/config_template.sh | tail -n1 >> ${installation_location}/new_config.sh
CPUs=$(nproc --all)
echo "procs=${CPUs}" >> ${installation_location}/new_config.sh
tail -n91 ${install_script_dir}/installation/config_template.sh | tail -n1 >> ${installation_location}/new_config.sh

# Copy all scripts from this fodler to install location
cp ${install_script_dir}/scripts/* ${installation_location}
chmod +X ${installation_location}/*


# Install databases
#Create database folder
if [[ ! -d ${databases} ]]; then
  echo "Creating ${databases} for databases"
  mkdir -p "${databases}"
fi
${installation_location}/database_checker.sh ${installation_location}/new_config.sh -i


#check if conda is installed already
conda_call_lines=$(conda list | wc -l)
if [[ "${conda_call_lines}" -gt 1 ]]; then
	:
else
	yes | ${install_script_dir}/installation/install_miniconda.sh
fi

home_dir=$($HOME)
echo "prefix: ${home_dir}/miniconda3/envs/py36_biopython_singularity" >> ${install_script_dir}/installation/py36_biopython_singularity.yml

conda create --name py36_biopython_singularity python=3.6 singularity biopython -y
