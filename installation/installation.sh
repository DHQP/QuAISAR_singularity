#! /usr/bin/env bash
#!/bin/bash -l


#
# Description: Script to prepare environment to be able to run quaisar pipeline
#
# Usage: ./installation.sh -i full_path_for_scripts -d full_path_to_download_databases_to -w full_path_of_where_to_store_output_of_runs
#
# Output location: Script, Database, and Output folders created and populated based on parameters
#
# Modules required: None
#
# v1.0.3 (03/29/2021)
#
# Created by Nick Vlachos (nvx4@cdc.gov)
#

#finds where script is, so it properly reference directories during install
install_script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd | rev | cut -d'/' -f2- | rev)"

echo "Current directory is ${install_script_dir}"

#  Function to print out help blurb
show_help () {
	echo "Usage: ./installation.sh -i full_path_for_scripts -d full_path_to_download_databases_to -w full_path_of_where_to_store_output_of_runs"
	echo " All paths need to be full, not relative"
}

# Parse command line options
options_found=0
while getopts ":h?i:d:w:" option; do
	options_found=$(( options_found + 1 ))
	case "${option}" in
		\?)
      echo "Invalid option found: ${OPTARG}"
      show_help
      exit 0
      ;;
		i)
      echo "Option -i triggered, argument = ${OPTARG}"
      installation_location=${OPTARG}
      ;;
		d)
      echo "Option -d triggered, argument = ${OPTARG}"
      databases=${OPTARG}
      ;;
		w)
      echo "Option -w triggered, argument = ${OPTARG}"
      working_directory=${OPTARG}
      ;;
		:)
      echo "Option -${OPTARG} requires as argument"
      ;;
		h)
      show_help
      exit 0
      ;;
	esac
done

# Show help info for when no options are given
if [[ "${options_found}" -eq 0 ]]; then
	echo "No options found"
	show_help
	exit 1
elif [ -z "${installation_location}" ]; then
	echo "Empty script location supplied to supplied to installation.sh, exiting"
	exit 3
elif [ -z "${databases}" ]; then
	echo "Empty database installation location supplied to installation.sh, exiting"
	exit 4
elif [ -z "${working_directory}" ]; then
  echo "Empty working directory for output supplied to installation.sh, exiting"
  exit 5
fi


echo -e "${installation_location}\n${databases}\n${working_directory}"

echo "Installing quaisar scripts to ${installation_location}"
if [[ ! -d ${installation_location} ]]; then
  echo "Creating ${installation_location}"
  mkdir -p ${installation_location}
fi

function test_go() {
  current_dir=$(pwd)
  echo "Testing Go!"
  echo -e "package main\n" > ${installation_location}/hello.go
  echo -e "import \"fmt\"\n" >> ${installation_location}/hello.go
  echo -e "func main() {\n\tfmt.Printf(\"hello, world\")\n}" >> ${installation_location}/hello.go
	cd ${installation_location}
	go build hello.go
	sleep 2
  go_test_output=$(./hello)
  if [[ "${go_test_output}" = "hello, world" ]]; then
    echo "Go is installed correctly, proceed!"
  else
    echo "Go is not installed, halting"
    exit 99
  fi
  cd ${current_dir}
}

function test_singularity() {
  sing_version=$(singularity --version)
  if [[ "${sing_version}" = "singularity version"* ]]; then
    v3plus=$(echo "${sing_version}" | cut -d' ' -f3 | cut -d'.' -f1)
    echo "${v3plus}"
    if [[ "${v3plus}" -ge 3 ]]; then
      echo "Singularity installed and ready to go"
    else
      echo "Your singularity version is too old, please upgrade to version 3 (or remove old versions and change parameters) before trying again"
      exit 100
    fi
  else
    echo "You do not have singularity installed, either choose a supported operating system, or preinstall before trying to install QuAISAR"
    exit 101
  fi
}

# function test_curl() {
# 	curled=$(curl -V)
# 	if [[ "${curled}" = "curl 7."* ]]; then
# 		echo "curl version 7.X installed, good to go"
# 	else
# 		echo "curl is NOT installed and is needed to download some databases, please install it using 'sudo apt get install curl' or 'sudo yum install curl'"
# 		exit 102
# 	fi
# }

# Do checks of singularity and Go to ensure they are installed already
test_go
test_singularity

echo "Go and Singularity are installed...proceed with installation"

echo "Checking for working directory of runs ${working_directory}"
if [[ ! -d ${working_directory} ]]; then
  echo "Creating ${working_directory}"
  mkdir -p ${working_directory}
fi

# Create primary config file
head -n25 ${install_script_dir}/installation/config_template.sh > ${installation_location}/new_config.sh
echo "output_dir=${working_directory}" >> ${installation_location}/new_config.sh
head -n27 ${install_script_dir}/installation/config_template.sh | tail -n1 >> ${installation_location}/new_config.sh
echo "src=${installation_location}" >> ${installation_location}/new_config.sh
head -n29 ${install_script_dir}/installation/config_template.sh | tail -n1 >> ${installation_location}/new_config.sh
echo "local_DBs=${databases}" >> ${installation_location}/new_config.sh
head -n31 ${install_script_dir}/installation/config_template.sh | tail -n1 >> ${installation_location}/new_config.sh
CPUs=$(nproc --all)
echo "procs=${CPUs}" >> ${installation_location}/new_config.sh
tail -n91 ${install_script_dir}/installation/config_template.sh >> ${installation_location}/new_config.sh
echo "source ${HOME}/miniconda3/etc/profile.d/conda.sh" >> ${installation_location}/new_config.sh

# Copy all scripts from this folder to install location
cp ${install_script_dir}/scripts/* ${installation_location}
rm ${installation_location}/config.sh
mv ${installation_location}/new_config.sh ${installation_location}/config.sh
echo "${installation_location}"
chmod +x ${installation_location}/*.sh ${installation_location}/*.py


# Install databases
#Create database folder
if [[ ! -d ${databases} ]]; then
  echo "Creating ${databases} for databases"
  mkdir -p "${databases}/BUSCO"
fi
echo "${install_script_dir}/scripts/database_checker.sh ${databases} -i"
${install_script_dir}/scripts/database_checker.sh ${databases} -i


#check if conda is installed already
conda_call_lines=$(conda list | wc -l)
if [[ "${conda_call_lines}" -gt 1 ]]; then
	:
else
	yes | "${install_script_dir}/installation/install_miniconda.sh"
fi

home_dir=$(echo $HOME)
echo "prefix: ${home_dir}/miniconda3/envs/py36_biopython" >> ${install_script_dir}/installation/py36_biopython.yml
echo 'export PATH=$PATH:${home_dir}/miniconda3/bin' >> ~/.bashrc
. "${home_dir}/.bashrc"

conda create --name py36_biopython python=3.6 biopython -y

echo -e "Installation complete. Before running the pipeline, open a new window or run the following:\n. ${home_dir}/.bashrc \n\n If you would like to view or edit parameters and settings for pipeline, please see ${installation_location}/config.sh"
