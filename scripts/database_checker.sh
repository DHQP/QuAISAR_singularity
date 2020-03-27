#!/bin/bash -l

#$ -o database_checker.out
#$ -e database_checker.err
#$ -N database_checker
#$ -cwd
#$ -q short.q

#Import the config file with shortcuts and settings
if [[ ! -f "./config.sh" ]]; then
	cp ./config_template.sh ./config.sh
fi
. ./config.sh

#
# Description: Script checks for all databases used by QuAISAR pipeline and sets up any missing ones
#
# Usage ./database_checker.sh -path_to_config_file [-i]
#	optional -i is for chencking AND installing, otherwise script just checks for existence..no downloading
# Output location: ${local_DBs} in config file
#
# Modules required: None
#
# v1.0 (03/04/2020)
#
# Created by Nick Vlachos (nvx4@cdc.gov)
#

if [[ ! -f "${1}" ]]; then
	echo "No config file...exiting"
else
	. "${1}"
fi

do_download="false"
# Checks for proper argumentation
if [[ $# -eq 0 ]]; then
	echo "No argument supplied to $0, exiting"
	exit 1
elif [[ -z "${1}" ]]; then
	echo "No config file...exiting"
	exit 1
elif [[ "${1}" = "-h" ]]; then
	echo "Usage ./database_checker.sh -path_to_config_file [-i]"
	echo "-i is too install databases, otherwise script just checks for existence"
	exit 0
elif [[ "${2}" == "-i" ]]; then
	do_download="true"
	#finds where script is it, so it properly reference directories during install
	current_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd | rev | cut -d'/' -f2- | rev)"
	if [[ ! -d ${current_dir}/installation ]]; then
		echo "Can not install databases. Try running ./database_checker.sh from the installation folder in the GIT repo"
		echo "${current_dir}"
	else
		do_download="true"
	fi
fi

# Shows where databases should be (installed)
echo "${local_DBs}"

missing_DBS=()

# Check for parent directory
if [[ ! -d ${local_DBs} ]]; then
	mkdir -p ${local_DBs}
fi

# Check for BUSCO
busco_taxa=(bacteria_odb10.2019-06-26 alteromonadales_odb10.2019-04-24 bacillales_odb10.2019-04-24 bacilli_odb10.2019-04-24 bacteroidetes_odb10.2019-04-24 betaproteobacteria_odb10.2019-04-24 burkholderiales_odb10.2019-04-24 campylobacterales_odb10.2019-04-24 clostridiales_odb10.2019-04-24 clostridia_odb10.2019-04-24 corynebacteriales_odb10.2019-04-24 enterobacterales_odb10.2019-04-24 epsilonproteobacteria_odb10.2019-04-24 firmicutes_odb10.2019-04-24 flavobacteriales_odb10.2019-04-24  flavobacteriia_odb10.2019-04-24 gammaproteobacteria_odb10.2019-04-24 lactobacillales_odb10.2019-04-24 neisseriales_odb10.2019-04-24 proteobacteria_odb10.2019-04-24 pseudomonadales_odb10.2019-04-24 xanthomonadales_odb10.2019-04-24 actinobacteria_class_odb10.2019-04-24 )

#echo "${#busco_taxa[@]}"

for odb_info in "${busco_taxa[@]}"; do
	# Check for top level bacteria database
	#echo ${odb_info}
	taxa=$(echo "$odb_info" | cut -d'_' -f1)
	db_date=$(echo "$odb_info" | cut -d'.' -f2)
	if [[ ! -d "${local_DBs}/BUSCO/${taxa}_odb10" ]]; then
		if [[ "${do_download}" = "true" ]]; then
			if [[ ! -d "${local_DBs}/BUSCO" ]]; then
				mkdir "${local_DBs}/BUSCO"
			fi
			cd "${local_DBs}/BUSCO"
			if [[ "${taxa}" == "actinobacteria" ]]; then
				taxa="actinobacteria_class"
			fi
			echo "Downloading latest BUSCO database for ${taxa} (wget http://busco-data.ezlab.org/v4/data/lineages/${taxa}_odb10.${db_date}.tar.gz)"
			wget "http://busco-data.ezlab.org/v4/data/lineages/${taxa}_odb10.${db_date}.tar.gz"
			# Dont know how to handle this one outlier (only one to specify a level in the filename) - ALl OUR bugs are in class Actinobacteria too
			if [[ "${taxa}" == "actinobacteria_class" ]]; then
				mv ${local_DBs}/BUSCO/actinobacteria_class_odb10.2019-04-24.tar.gz ${local_DBs}/BUSCO/actinobacteria_odb10.2019-04-24.tar.gz
			fi
		else
			echo "Missing latest BUSCO database for ${taxa}"
			missing_DBS=("${missing_DBS[@]}" "BUSCO-${taxa}")
		fi
	else
		echo "BUSCO has latest ${taxa}_odb10 as of 3/15/2020"
	fi
done
find ${local_DBs}/BUSCO/ -name '*.gz' -exec tar xzf {} \;
find ${local_DBs}/BUSCO/ -name '*.gz' -exec rm {} \;

##### Currently down.....and has been a while
# Check to see if gottcha database is installed
if [[ ! -d "${local_DBs}/gottcha" ]]; then
	if [[ "${do_download}" = "true" ]]; then
		cd "${local_DBs}"
		# Original LANL hosted address that has been down a good while
	 	#wget -P "${local_DBs}/gottcha" "https://edge-dl.lanl.gov/gottcha/GOTTCHA_database_v20150825/GOTTCHA_BACTERIA_c4937_k24_u30_xHUMAN3x.species.tar.gz"
		# Temporary mirror until original is fixed
		echo "Downloading latest gottcha database (wget https://zenodo.org/record/819341/files/gottcha_bac_arc_v1.tar.gz)"
		wget "https://zenodo.org/record/819341/files/gottcha_bac_arc_v1.tar.gz"
		tar xzf gottcha_bac_arc_v1.tar.gz
		rm gottcha_bac_arc_v1.tar.gz
		mv gottcha/gottcha_db ./
		rm -r gottcha
		mv gottcha_db gottcha
		rm gottcha.dbprofile.out
		# Need to find sa place to host genus_Lookup.tar.gz
	else
		echo "Missing gottcha database"
		missing_DBS=("${missing_DBS[@]}" "gottcha")
	fi
else
	echo "gottcha database installed"
fi

# Check to see if kraken mini database is installed
if [[ ! -d "${local_DBs}/kraken" ]]; then
	if [[ "${do_download}" = "true" ]]; then
		mkdir "${local_DBs}/kraken"
		cd "${local_DBs}/kraken"
		echo "Downloading latest (mini)kraken database (wget https://ccb.jhu.edu/software/kraken/dl/minikraken_20171019_4GB.tgz)"
		wget "https://ccb.jhu.edu/software/kraken/dl/minikraken_20171019_4GB.tgz"
		tar xzf minikraken_20171019_4GB.tgz
		rm minikraken_20171019_4GB.tgz
	else
		echo "Missing latest kraken database"
		missing_DBS=("${missing_DBS[@]}" "kraken")
	fi
else
	echo "kraken database is installed"
fi

# All other databases will need to be hosted somehwere before being able to be checked/updated. Currently they are included in the Docker image

# ANI sketch file / aniDB (150MBs/1.3Gbs)
if [[ ! -d "${local_DBs}/ANI" ]]; then
	#cp -r /container_DBs/ANI ${local_DBs}
	if [[ "${do_download}" = "true" ]]; then
		mkdir "${local_DBs}/ANI"
		echo "Copying latest REFSEQ sketch database (ANI)"
		cp ${current_dir}/databases/ANI ${local_DBs}
	else
		echo "Missing latest REFSEQ sketch database (ANI)"
		missing_DBS=("${missing_DBS[@]}" "REFSEQ-ANI")
	fi
else
	echo "ANI REFSEQ sketch database installed"
fi

# star (6 Mbs)
if [[ ! -d "${local_DBs}/star" ]]; then
	#cp -r /container_DBs/star ${local_DBs}
	if [[ "${do_download}" = "true" ]]; then
		echo "Copying latest NAR-AR database"
		cp -r ${current_dir}/included_databases/star ${local_DBs}
	else
		echo "Missing latest NAR-AR database"
		missing_DBS=("${missing_DBS[@]}" "NAR-AR")
	fi
else
	echo "NAR-AR database installed"
fi


singularities=(bbtools.simg blast-2.9.0-docker.img bowtie2-2.2.9-biocontainers.simg cSSTAR.simg entrez_taxon.simg GAMA_quaisar.simg gottcha.simg plasmidFinder_with_DB.simg QUAST5.simg srst2.simg)

for simage in "${singularities[@]}"; do
	# custom singularity images (3.6GBs)
	if [[ ! -f "${local_DBs}/singularities/${simage}" ]]; then
		#cp -r /container_DBs/custom_singularities ${local_DBs}
		if [[ "${do_download}" = "true" ]]; then
			if [[ ! -d "${local_DBs}/singularities" ]]; then
				mkdir "${local_DBs}/singularities"
			fi
			echo "Copying custom singularity image ${simage}"
			cp ${current_dir}/included_databases/singularities/${simage} ${local_DBs}/singularities
		else
			echo "Missing custom singularity image ${simage}"
			missing_DBS=("${missing_DBS[@]}" "singularities-${simage}")
		fi
	else
		echo "custom singularity image ${simage} installed"
	fi
done

if [[ ! -f "${local_DBs}/MMB_Bugs.txt" ]]; then
	#cp -r /container_DBs/MMB_Bugs.txt ${local_DBs}
	if [[ "${do_download}" = "true" ]]; then
		echo "Copying MMB_bugs"
		cp ${current_dir}/included_databases/MMB_Bugs.txt ${local_DBs}
	else
		echo "Missing MMB_Bugs"
		missing_DBS=("${missing_DBS[@]}" "MMB_Bugs")
	fi
else
	echo "MMB_Bugs installed"
fi

if [[ ! -f "${local_DBs}/taxes.csv" ]]; then
	#cp -r /container_DBs/taxes.csv ${local_DBs}
	if [[ "${do_download}" = "true" ]]; then
		echo "Copying taxes"
		cp ${current_dir}/included_databases/taxes.csv ${local_DBs}
	else
		echo "Missing taxes"
		missing_DBS=("${missing_DBS[@]}" "taxes")
	fi
else
	echo "taxes installed"
fi

if [[ ! -f "${local_DBs}/phiX.fasta" ]]; then
	#cp -r /container_DBs/phiX.fasta ${local_DBs}
	if [[ "${do_download}" = "true" ]]; then
		echo "Copying phiX.fasta"
		cp ${current_dir}/included_databases/phiX.fasta ${local_DBs}
	else
		echo "Missing phiX"
		missing_DBS=("${missing_DBS[@]}" "phiX")
	fi
else
	echo "phiX installed"
fi

if [[ ! -f "${local_DBs}/adapters.fasta" ]]; then
	#cp -r /container_DBs/adapters.fasta ${local_DBs}
	if [[ "${do_download}" = "true" ]]; then
		echo "Copying adapters.fasta"
		cp ${current_dir}/included_databases/adapters.fasta ${local_DBs}
	else
		echo "Missing adapters"
		missing_DBS=("${missing_DBS[@]}" "adapters")
	fi
else
	echo "adapters installed"
fi

if [[ ! -d "${local_DBs}/pubmlsts" ]]; then
	#cp -r /container_DBs/pubmlsts ${local_DBs}
	if [[ "${do_download}" = "true" ]]; then
		echo "Copying pubMLST"
		cp ${current_dir}/included_databases/pubmlsts ${local_DBs}
	else
		echo "Missing pubMLST"
		missing_DBS=("${missing_DBS[@]}" "pubMLST")
	fi
else
	echo "pubMLST installed"
fi

echo "There are ${#missing_DBS[@]} missing databases (${missing_DBS[@]})"
