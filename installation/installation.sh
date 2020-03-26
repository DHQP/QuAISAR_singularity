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


# need to add to bottom of yaml with proper home location
installation_location=${1}
databases=${2}
working_directory=${3}

echo "Installing quaisar scripts to ${installation_location}"
if [[ ! -d ${installation_location} ]]; then
  echo "Creating ${installation}"
  mkdir -p ${installation_location}
fi

# Create primary config file
head -n25 ./config_template.sh > ${installation_location}/config.sh
echo "output_dir=${working_directory}" >> ${installation_location}/config.sh
head -n27 ./config_template.sh | tail -n1 >> ${installation_location}/config.sh
echo "src=${installation_location}" >> ${installation_location}/config.sh
head -n29 ./config_template.sh | tail -n1 >> ${installation_location}/config.sh
echo "local_DBs=${databases}" >> ${installation_location}/config.sh
head -n31 ./config_template.sh | tail -n1 >> ${installation_location}/config.sh
CPUs=$(nproc --all)
echo "procs=${CPUs}" >> ${installation_location}/config.sh
tail -n91 ./config_template.sh | tail -n1 >> ${installation_location}/config.sh

# Install databases
${source}/database_checker.sh $config_file -i

# Copy all scripts from this fodler to install location
cp ../scripts/* ${installation_location}


#check if conda is installed already
conda_call_lines=$(conda list | wc -l)
if [[ "${conda_call_lines}" -gt 1 ]]; then
	:
else
	yes | ${source}/installation/install_miniconda.sh
fi

home_dir=$($HOME)
echo "prefix: ${home_dir}/miniconda3/envs/py36_biopython_singularity" >> ./py36_biopython_singularity.yml

conda create --name py36_biopython_singularity python=3.6 singularity biopython -y
