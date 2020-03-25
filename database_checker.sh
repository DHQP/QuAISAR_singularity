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
fi

# Shows where databases should be (installed)
echo "${local_DBs}"

missing_DBS=()

# Check for parent directory
if [[ ! -d ${local_DBs} ]]; then
	mkdir -p ${local_DBs}
fi

# Check for BUSCO
busco_taxa=(bacteria_odb10.2019-06-26 alteromonadales_odb10.2019-04-24 bacillales_odb10.2019-04-24 bacilli_odb10.2019-04-24 bacteroidetes_odb10.2019-04-24 betaproteobacteria_odb10.2019-04-24 burkholderiales_odb10.2019-04-24 campylobacterales_odb10.2019-04-24 clostridiales_odb10.2019-04-24 clostridia_odb10.2019-04-24 corynebacteriales_odb10.2019-04-24 enterobacterales_odb10.2019-04-24 epsilonproteobacteria_odb10.2019-04-24 firmicutes_odb10.2019-04-24 flavobacteriales_odb10.2019-04-24 flavobacteriales_odb10.2019-04-24 flavobacteriia_odb10.2019-04-24 gammaproteobacteria_odb10.2019-04-24 lactobacillales_odb10.2019-04-24 neisseriales_odb10.2019-04-24 proteobacteria_odb10.2019-04-24 pseudomonadales_odb10.2019-04-24 xanthomonadales_odb10.2019-04-24 actinobacteria_class_odb10.2019-04-24)

echo "${#busco_taxa[@]}"

for odb_info in "${busco_taxa[@]}"; do
	# Check for top level bacteria database
	echo ${odb_info}
	taxa=$(echo "$odb_info" | cut -d'_' -f1)
	db_date=$(echo "$odb_info" | cut -d'.' -f2)
	if [[ ! -d "${local_DBs}/${taxa}_odb10" ]]; then
		if [[ "${do_download}" = "true" ]]; then
			if [[ ! -d "${local_DBs}/BUSCO" ]]; then
				mkdir "${local_DBs}/BUSCO"
				cd "${local_DBs}/BUSCO"
			fi
			echo "Downloading latest BUSCO database for ${taxa} (wget http://busco-data.ezlab.org/v4/data/lineages/${taxa}_odb10.${db_date}.tar.gz)"
			wget "http://busco-data.ezlab.org/v4/data/lineages/${taxa}_odb10.${db_date}.tar.gz"
		else
			echo "Missing latest BUSCO database for ${taxa}"
			missing_DBS=("${missing_DBS[@]}" "BUSO-${taxa}")
		fi
	else
		echo "BUSCO has latest ${taxa}_odb10 as of 3/15/2020"
	fi
done

	# # Check for top level bacteria database
	# if [[ ! -d "${local_DBs}/bacteria_odb10" ]]; then
	# 	wget "http://busco-data.ezlab.org/v4/data/lineages/bacteria_odb10.2019-06-26.tar.gz"
	# else
	# 	echo "BUSCO has latest Bacteria as of 3/15/2020"
	# fi
	# # Check for Order level Alteromonadales database
	# if [[ ! -d "${local_DBs}/alteromonadales_odb10" ]]; then
	# 	wget "http://busco-data.ezlab.org/v4/data/lineages/alteromonadales_odb10.2019-04-24.tar.gz"
	# fi
	# # Check for Order level Bacillales database
	# if [[ ! -d "${local_DBs}/bacillales_odb10" ]]; then
	# 	wget "http://busco-data.ezlab.org/v4/data/lineages/bacillales_odb10.2019-04-24.tar.gz"
	# fi
	# # Check for Class level Bacilli database
	# if [[ ! -d "${local_DBs}/bacilli_odb10" ]]; then
	# 	wget "http://busco-data.ezlab.org/v4/data/lineages/bacilli_odb10.2019-04-24.tar.gz"
	# fi
	# # Check for Phylum level Bacteroidetes database
	# if [[ ! -d "${local_DBs}/bacteroidetes_odb10" ]]; then
	# 	wget "http://busco-data.ezlab.org/v4/data/lineages/bacteroidetes_odb10.2019-04-24.tar.gz"
	# fi
	# # Check for Class level Betaproteobacteria database
	# if [[ ! -d "${local_DBs}/betaproteobacteria_odb10" ]]; then
	# 	wget "http://busco-data.ezlab.org/v4/data/lineages/betaproteobacteria_odb10.2019-04-24.tar.gz"
	# fi
	# # Check for Order level Burkholderiales database
	# if [[ ! -d "${local_DBs}/burkholderiales_odb10" ]]; then
	# 	wget "http://busco-data.ezlab.org/v4/data/lineages/burkholderiales_odb10.2019-04-24.tar.gz"
	# fi
	# # Check for Order level Campylobacterales database
	# if [[ ! -d "${local_DBs}/campylobacterales_odb10" ]]; then
	# 	wget "http://busco-data.ezlab.org/v4/data/lineages/campylobacterales_odb10.2019-04-24.tar.gz"
	# fi
	# # Check for Order level Clostridiales database
	# if [[ ! -d "${local_DBs}/clostridiales_odb10" ]]; then
	# 	wget "http://busco-data.ezlab.org/v4/data/lineages/clostridiales_odb10.2019-04-24.tar.gz"
	# fi
	# # Check for Class level Clostridia database
	# if [[ ! -d "${local_DBs}/clostridia_odb10" ]]; then
	# 	wget "http://busco-data.ezlab.org/v4/data/lineages/clostridia_odb10.2019-04-24.tar.gz"
	# fi
	# # Check for Order level Corynebacteriales database
	# if [[ ! -d "${local_DBs}/corynebacteriales_odb10" ]]; then
	# 	wget "http://busco-data.ezlab.org/v4/data/lineages/corynebacteriales_odb10.2019-04-24.tar.gz"
	# fi
	# # Check for Order level Enterobacterales database
	# if [[ ! -d "${local_DBs}/enterobacterales_odb10" ]]; then
	# 	wget "http://busco-data.ezlab.org/v4/data/lineages/enterobacterales_odb10.2019-04-24.tar.gz"
	# fi
	# # Check for Class level Epsilonproteobacteria
	# if [[ ! -d "${local_DBs}/epsilonproteobacteria_odb10" ]]; then
	# 	wget "http://busco-data.ezlab.org/v4/data/lineages/epsilonproteobacteria_odb10.2019-04-24.tar.gz"
	# fi
	# # Check for Phylum level Firmicutes database
	# if [[ ! -d "${local_DBs}/firmicutes_odb10" ]]; then
	# 	wget "http://busco-data.ezlab.org/v4/data/lineages/firmicutes_odb10.2019-04-24.tar.gz"
	# fi
	# # Check for Order level Flavobacteriales database
	# if [[ ! -d "${local_DBs}/flavobacteriales_odb10" ]]; then
	# 	wget "http://busco-data.ezlab.org/v4/data/lineages/flavobacteriales_odb10.2019-04-24.tar.gz"
	# fi
	# # Check for Class level Flavobacteriia database
	# if [[ ! -d "${local_DBs}/flavobacteriia_odb10" ]]; then
	# 	wget "http://busco-data.ezlab.org/v4/data/lineages/flavobacteriia_odb10.2019-04-24.tar.gz"
	# fi
	# # Check for Class level Gammaproteobacteria database
	# if [[ ! -d "${local_DBs}/gammaproteobacteria_odb10" ]]; then
	# 	wget "http://busco-data.ezlab.org/v4/data/lineages/gammaproteobacteria_odb10.2019-04-24.tar.gz"
	# fi
	# # Check for Order level lactobacillales database
	# if [[ ! -d "${local_DBs}/lactobacillales_odb10" ]]; then
	# 	wget "http://busco-data.ezlab.org/v4/data/lineages/lactobacillales_odb10.2019-04-24.tar.gz"
	# fi
	# # Check for Order level Neisseriales database
	# if [[ ! -d "${local_DBs}/neisseriales_odb10" ]]; then
	# 	wget "http://busco-data.ezlab.org/v4/data/lineages/neisseriales_odb10.2019-04-24.tar.gz"
	# fi
	# # Check for Phylum level Proteobacteria database
	# if [[ ! -d "${local_DBs}/proteobacteria_odb10" ]]; then
	# 	wget "http://busco-data.ezlab.org/v4/data/lineages/proteobacteria_odb10.2019-04-24.tar.gz"
	# fi
	# # Check for Order level Pseudomonadales database
	# if [[ ! -d "${local_DBs}/pseudomonadales_odb10" ]]; then
	# 	wget "http://busco-data.ezlab.org/v4/data/lineages/pseudomonadales_odb10.2019-04-24.tar.gz"
	# fi
	# # Check for Order level Xanthomonadales database
	# if [[ ! -d "${local_DBs}/xanthomonadales_odb10" ]]; then
	# 	wget "http://busco-data.ezlab.org/v4/data/lineages/xanthomonadales_odb10.2019-04-24.tar.gz"
	# fi
	# # Check for Class level Actinobacteria database. This is the only conflict database as there is also a phylum level Actinobacteria which is not going to be used.
	# if [[ ! -d "${local_DBs}/actinobacteria_class_odb10" ]]; then
	# 	wget "http://busco-data.ezlab.org/v4/data/lineages/actinobacteria_class_odb10.2019-04-24.tar.gz"
	# fi
	for file in ${local_DBs}/BUSCO/*.gz; do
		tar xzf ${file}
	done
	rm *tar.gz
fi

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
		echo "Copying latest REFSEQ sketch database (ANI)"
		cp ${src}/databases/ANI ${local_DBs}
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
		cp ${src}/databases/star ${local_DBs}
	else
		echo "Missing latest NAR-AR database"
		missing_DBS=("${missing_DBS[@]}" "NAR-AR")
	fi
else
	echo "NAR-AR database installed"
fi


singularities=(bbtools blast-2.9.0-docker bowtie2-2.2.9-biocontainers cSSTAR entrez_taxon GAMA_quaisar gottcha plasmidFinder_with_DB QUAST5 srst2)

for simage in "${singularities[@]}"; do
	# custom singularity images (3.6GBs)
	if [[ ! -f "${local_DBs}/custom_singularities/${simage}.simg" ]]; then
		#cp -r /container_DBs/custom_singularities ${local_DBs}
		if [[ "${do_download}" = "true" ]]; then
			if [[ ! -d "${local_DBs}/custom_singularities" ]]; then
				mkdir "${local_DBs}/custom_singularities"
			fi
			echo "Copying custom singularity image ${simage}.simg"
			cp ${src}/databases/singularities/${simage}.simg ${local_DBs}/custom_singularities
		else
			echo "Missing custom singularity image ${simage}.simg"
			missing_DBS=("${missing_DBS[@]}" "singularities-${simage}")
		fi
	else
		echo "custom singularity image ${simage}.simg installed"
	fi
done

if [[ ! -f "${local_DBs}/MMB_Bugs.txt" ]]; then
	#cp -r /container_DBs/MMB_Bugs.txt ${local_DBs}
	if [[ "${do_download}" = "true" ]]; then
		echo "Copying MMB_bugs"
		cp ${src}/databases/MMB_bugs.txt ${local_DBs}
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
		cp ${src}/databases/taxes.csv ${local_DBs}
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
		cp ${src}/databases/phiX.fasta ${local_DBs}
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
		cp ${src}/databases/adapters.fasta ${local_DBs}
	else
		echo "Missing adapters"
		missing_DBS=("${missing_DBS[@]}" "adapters")
	fi
else
	echo "adapters installed"
fi

if [[ ! -d "${local_DBs}/pubMLSTs" ]]; then
	#cp -r /container_DBs/pubMLSTs ${local_DBs}
	if [[ "${do_download}" = "true" ]]; then
		echo "Copying pubMLST"
		cp ${src}/databases/pubMLST ${local_DBs}
	else
		echo "Missing pubMLST"
		missing_DBS=("${missing_DBS[@]}" "pubMLST")
	fi
else
	echo "pubMLST installed"
fi
