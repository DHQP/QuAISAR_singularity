#! /usr/bin/env bash
#!/bin/sh -l


#
# Description: Script to prepare environment to be able to run quaisar pipeline
#
# Usage: ./quaisar_installation.sh OS_Type location_for_scriopts location_to_download_databases_to location_of_where_to_store_output
#
# Output location: No output created
#
# Modules required: None
#
# v1.0.1 (04/01/2020)
#
# Created by Nick Vlachos (nvx4@cdc.gov)
#

#finds where script is it, so it properly reference directories during install
install_script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd | rev | cut -d'/' -f2- | rev)"

echo "${install_script_dir}"

# Checking for proper number of arguments from command line
if [[ $# -lt 1  || $# -gt 4 ]]; then
  echo "Usage: ./installation.sh OS-Type script_installation_location database_installation_location Working_directory_of_output_from_runs"
	echo "You have used $# args"
  echo "Options for OS-type are 0-Other/No Singularity Install needed 1-Ubuntu/Debian Family, 2-CentOS/Red Hat Family"
  exit 3
fi

echo "Installing quaisar scripts to ${installation_location}"
if [[ ! -d ${installation_location} ]]; then
  echo "Creating ${installation_location}"
  mkdir -p ${installation_location}
fi

# Checks for proper argumentation
if [[ $# -eq 0 ]]; then
	echo "No argument supplied to $0, exiting"
	exit 1
elif [[ "$1" = "-h" ]]; then
	echo "Usage: ./installation.sh script_installation_location database_installation_location Working_directory_of_output_from_runs"
  exit 0
elif [ -z "$1" ]; then
  echo "Empty OS_type supplied to supplied to $0, exiting"
  exit 1
elif [ -z "$2" ]; then
	echo "Empty script location supplied to supplied to $0, exiting"
	exit 1
elif [ -z "$3" ]; then
	echo "Empty database installation location supplied to $0, exiting"
	exit 1
elif [ -z "$4" ]; then
  echo "Empty working directory for output supplied to $0, exiting"
  exit 1
fi

function test_go() {
  current_dir=$(pwd)
  "Testing Go!"
  echo -e "package main\n" > ${installation_location}/hello.go
  echo -e "import fmt\n" >> ${installation_location}/hello.go
  echo -e "func main() {\nfmt.Printf(\"hello, world\")\n}" >> ${installation_location}/hello.go
  go build ${installation_location}/hello.go
  cd ${installation_location}
  go_test_output=$(go hello)
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
    if [[ "${v3plus}" -gt 3 ]]; then
      echo "Singularity installed and ready to go"
    else
      echo "Your singularity version is too old, please upgrade to version 3 (or remove old versions and change parameters) before trying again"
      exit
    fi
  else
    echo "You do not have singularity installed, either choose a supported operating system, or preinstall before trying to install QuAISAR"
    exit 101
  fi
}

# need to add to bottom of yaml with proper home location
OS_type=${1}
installation_location=${2}
databases=${3}
working_directory=${4}

#Check if Go is installed already






echo "Checking for working directory of runs ${working_directory}"
if [[ ! -d ${working_directory} ]]; then
  echo "Creating ${working_directory}"
  mkdir -p ${working_directory}
fi

if [[ ${OS_type} -eq 1 ]]; then
  echo "Installing pre-dependencies"
  sudo apt-get update && sudo apt-get install -y \
    build-essential \
    libssl-dev \
    uuid-dev \
    libgpgme11-dev \
    squashfs-tools \
    libseccomp-dev \
    wget \
    pkg-config \
    git \
    cryptsetup
  echo "Installing Go"
  curl -O https://storage.googleapis.com/golang/go1.14.2.linux-amd64.tar.gz --output ${installation_location}
  tar -C /usr/local -xvf go1.14.2.linux-amd64.tar.gz
  echo "export PATH=$PATH:/usr/local/go/bin" >> /etc/profile
  test_go
  echo "Installing Singularity"
  export VERSION=3.5.2 && # adjust this as necessary \
    wget https://github.com/sylabs/singularity/releases/download/v${VERSION}/singularity-${VERSION}.tar.gz && \
    tar -C "${installation_location}" -xzf singularity-${VERSION}.tar.gz && \
    cd ${installation_location}/singularity
  ./mconfig && \
    make -c builddir && \
    sudo make -C builddir install.
  test_singularity

#Install Singularity on Redhat based systems
elif [[ ${OS_type} -eq 2]]; then

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

# Copy all scripts from this fodler to install location
cp ${install_script_dir}/scripts/* ${installation_location}
rm ${installation_location}/config.sh
mv ${installation_location}/new_config.sh ${installation_location}/config.sh
chmod +x ${installation_location}/*

echo "Exiting prematurely to prevent constant installations"

# Install databases
#Create database folder
if [[ ! -d ${databases} ]]; then
  echo "Creating ${databases} for databases"
  mkdir -p "${databases}/BUSCO"
fi
${install_script_dir}/scripts/database_checker.sh ${installation_location}/config.sh -i


#check if conda is installed already
conda_call_lines=$(conda list | wc -l)
if [[ "${conda_call_lines}" -gt 1 ]]; then
	:
else
	yes | ${install_script_dir}/installation/install_miniconda.sh
fi

home_dir=$(echo $HOME)
echo "prefix: ${home_dir}/miniconda3/envs/py36_biopython" >> ${install_script_dir}/installation/py36_biopython.yml

. "${home_dir}/.bashrc"

conda create --name py36_biopython python=3.6 biopython -y

echo -e "Installation complete. Before running the pipeline, open a new window or run the following:\n. ${home_dir}/.bashrc \n\n If you would like to view or edit parameters and settings for pipeline, please see ${installation_location}/config.sh"
