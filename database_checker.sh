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
# Usage ./database_checker.sh -path_to_config_file
#
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

echo "${local_DBs}"

exit


# Check for parent directory
if [[ ! -d ${local_DBs} ]]; then
	mkdir -p ${local_DBs}
fi

# Check for BUSCO
if [[ ! -d "${local_DBs}/BUSCO" ]]; then
	mkdir -p "${local_DBs}/BUSCO"
	# Check for top level bacteria database
	if [[ ! -d "${local_DBs}/bacteria_odb10.2019-06-26" ]]; then
		wget -P "${local_DBs}/BUSCO" "http://busco-data.ezlab.org/v4/data/lineages/bacteria_odb10.2019-06-26.tar.gz"
	fi
	# Check for Order level Alteromonadales database
	if [[ ! -d "${local_DBs}/bacteria_odb10.2019-06-26" ]]; then
		wget -P "${local_DBs}/BUSCO" "http://busco-data.ezlab.org/v4/data/lineages/alteromonadales_odb10.2019-04-24.tar.gz"
	fi
	# Check for Order level Bacillales database
	if [[ ! -d "${local_DBs}/bacteria_odb10.2019-06-26" ]]; then
		wget -P "${local_DBs}/BUSCO" "http://busco-data.ezlab.org/v4/data/lineages/bacillales_odb10.2019-04-24.tar.gz"
	fi
	# Check for Class level Bacilli database
	if [[ ! -d "${local_DBs}/bacteria_odb10.2019-06-26" ]]; then
		wget -P "${local_DBs}/BUSCO" "http://busco-data.ezlab.org/v4/data/lineages/bacilli_odb10.2019-04-24.tar.gz"
	fi
	# Check for Phylum level Bacteroidetes database
	if [[ ! -d "${local_DBs}/bacteria_odb10.2019-06-26" ]]; then
		wget -P "${local_DBs}/BUSCO" "http://busco-data.ezlab.org/v4/data/lineages/bacteroidetes_odb10.2019-04-24.tar.gz"
	fi
	# Check for Class level Betaproteobacteria database
	if [[ ! -d "${local_DBs}/bacteria_odb10.2019-06-26" ]]; then
		wget -P "${local_DBs}/BUSCO" "http://busco-data.ezlab.org/v4/data/lineages/betaproteobacteria_odb10.2019-04-24.tar.gz"
	fi
	# Check for Order level Burkholderiales database
	if [[ ! -d "${local_DBs}/bacteria_odb10.2019-06-26" ]]; then
		wget -P "${local_DBs}/BUSCO" "http://busco-data.ezlab.org/v4/data/lineages/burkholderiales_odb10.2019-04-24.tar.gz"
	fi
	# Check for Order level Campylobacterales database
	if [[ ! -d "${local_DBs}/bacteria_odb10.2019-06-26" ]]; then
		wget -P "${local_DBs}/BUSCO" "http://busco-data.ezlab.org/v4/data/lineages/campylobacterales_odb10.2019-04-24.tar.gz"
	fi
	# Check for Order level Clostridiales database
	if [[ ! -d "${local_DBs}/bacteria_odb10.2019-06-26" ]]; then
		wget -P "${local_DBs}/BUSCO" "http://busco-data.ezlab.org/v4/data/lineages/clostridiales_odb10.2019-04-24.tar.gz"
	fi
	# Check for Class level Clostridia database
	if [[ ! -d "${local_DBs}/bacteria_odb10.2019-06-26" ]]; then
		wget -P "${local_DBs}/BUSCO" "http://busco-data.ezlab.org/v4/data/lineages/clostridia_odb10.2019-04-24.tar.gz"
	fi
	# Check for Order level Corynebacteriales database
	if [[ ! -d "${local_DBs}/bacteria_odb10.2019-06-26" ]]; then
		wget -P "${local_DBs}/BUSCO" "http://busco-data.ezlab.org/v4/data/lineages/corynebacteriales_odb10.2019-04-24.tar.gz"
	fi
	# Check for Order level Enterobacterales database
	if [[ ! -d "${local_DBs}/bacteria_odb10.2019-06-26" ]]; then
		wget -P "${local_DBs}/BUSCO" "http://busco-data.ezlab.org/v4/data/lineages/enterobacterales_odb10.2019-04-24.tar.gz"
	fi
	# Check for Class level Epsilonproteobacteria
	if [[ ! -d "${local_DBs}/bacteria_odb10.2019-06-26" ]]; then
		wget -P "${local_DBs}/BUSCO" "http://busco-data.ezlab.org/v4/data/lineages/epsilonproteobacteria_odb10.2019-04-24.tar.gz"
	fi
	# Check for Phylum level Firmicutes database
	if [[ ! -d "${local_DBs}/firmicutes_odb10.2019-04-24" ]]; then
		wget -P "${local_DBs}/BUSCO" "http://busco-data.ezlab.org/v4/data/lineages/firmicutes_odb10.2019-04-24.tar.gz"
	fi
	# Check for Order level Flavobacteriales database
	if [[ ! -d "${local_DBs}/flavobacteriales_odb10.2019-04-24" ]]; then
		wget -P "${local_DBs}/BUSCO" "http://busco-data.ezlab.org/v4/data/lineages/flavobacteriales_odb10.2019-04-24.tar.gz"
	fi
	# Check for Class level Flavobacteriia database
	if [[ ! -d "${local_DBs}/flavobacteriia_odb10.2019-04-24" ]]; then
		wget -P "${local_DBs}/BUSCO" "http://busco-data.ezlab.org/v4/data/lineages/flavobacteriia_odb10.2019-04-24.tar.gz"
	fi
	# Check for Class level Gammaproteobacteria database
	if [[ ! -d "${local_DBs}/gammaproteobacteria_odb10.2019-04-24" ]]; then
		wget -P "${local_DBs}/BUSCO" "http://busco-data.ezlab.org/v4/data/lineages/gammaproteobacteria_odb10.2019-04-24.tar.gz"
	fi
	# Check for Order level lactobacillales database
	if [[ ! -d "${local_DBs}/lactobacillales_odb10.2019-04-24" ]]; then
		wget -P "${local_DBs}/BUSCO" "http://busco-data.ezlab.org/v4/data/lineages/lactobacillales_odb10.2019-04-24.tar.gz"
	fi
	# Check for Order level Neisseriales database
	if [[ ! -d "${local_DBs}/neisseriales_odb10.2019-04-24" ]]; then
		wget -P "${local_DBs}/BUSCO" "http://busco-data.ezlab.org/v4/data/lineages/neisseriales_odb10.2019-04-24.tar.gz"
	fi
	# Check for Phylum level Proteobacteria database
	if [[ ! -d "${local_DBs}/proteobacteria_odb10.2019-04-24" ]]; then
		wget -P "${local_DBs}/BUSCO" "http://busco-data.ezlab.org/v4/data/lineages/proteobacteria_odb10.2019-04-24.tar.gz"
	fi
	# Check for Order level Pseudomonadales database
	if [[ ! -d "${local_DBs}/pseudomonadales_odb10.2019-04-24" ]]; then
		wget -P "${local_DBs}/BUSCO" "http://busco-data.ezlab.org/v4/data/lineages/pseudomonadales_odb10.2019-04-24.tar.gz"
	fi
	# Check for Order level Xanthomonadales database
	if [[ ! -d "${local_DBs}/xanthomonadales_odb10.2019-04-24" ]]; then
		wget -P "${local_DBs}/BUSCO" "http://busco-data.ezlab.org/v4/data/lineages/xanthomonadales_odb10.2019-04-24.tar.gz"
	fi
	# Check for Class level Actinobacteria database. This is the only conflict database as there is also a phylum level Actinobacteria which is not going to be used.
	if [[ ! -d "${local_DBs}/actinobacteria_class_odb10.2019-04-24" ]]; then
		wget -P "${local_DBs}/BUSCO" "http://busco-data.ezlab.org/v4/data/lineages/actinobacteria_class_odb10.2019-04-24.tar.gz"
	fi
	chmod 777 ${local_DBs}/BUSCO/*.tar.gz
	tar -vzfx ${local_DBs}/BUSCO/*.tar.gz
fi

##### Currently down.....and has been a while
# Check to see if gottcha database is installed
 if [[ ! -d "${local_DBs}/gottcha" ]]; then
 	mkdir -p "${local_DBs}/gottcha"
	# Original LANL hosted address that has been down a good while
 	#wget -P "${local_DBs}/gottcha" "https://edge-dl.lanl.gov/gottcha/GOTTCHA_database_v20150825/GOTTCHA_BACTERIA_c4937_k24_u30_xHUMAN3x.species.tar.gz"
	# Temporary mirror until original is fixed
	wget -O ${local_DBs}/gottcha/GOTTCHA_BACTERIA_c4937_k24_u30_xHUMAN3x.species.tar.gz "${local_DBs}/gottcha" "https://zenodo.org/record/819341/files/gottcha_bac_arc_v1.tar.gz?download=1"
	chmod 777 "${local_DBs}/gottcha/GOTTCHA_BACTERIA_c4937_k24_u30_xHUMAN3x.species.tar.gz"
	tar -xvzf "${local_DBs}/gottcha/GOTTCHA_BACTERIA_c4937_k24_u30_xHUMAN3x.species.tar.gz"
 fi

# # Possible solution to gottchaV1 not working
# if [[ ! -d "${local_DBs}/gottcha2" ]]; then
# 	mkdir "${local_DBs}/gottcha2"
# 	wget -O ${local_DBs}/gottcha2/taxdump.tar.gz https://edge-dl.lanl.gov/GOTTCHA2/RefSeq-Release90/taxdump.tar.gz
# 	wget https://edge-dl.lanl.gov/GOTTCHA2/RefSeq-Release90/RefSeq-r90.cg.BacteriaArchaeaViruses.species.fna.tar
# 	tar -xf RefSeq-r90.cg.BacteriaArchaeaViruses.species.fna.tar -C ${local_DBs}/gottcha2
# fi

# Check to see if kraken mini database is installed
if [[ ! -d "${local_DBs}/kraken" ]]; then
	mkdir -p "${local_DBs}/kraken"
	wget -P "${local_DBs}/kraken" "https://ccb.jhu.edu/software/kraken/dl/minikraken_20171019_4GB.tgz"
	chmod 777 "${local_DBs}/kraken/minikraken_20171019_4GB.tgz"
	tar -xzvf "${local_DBs}/kraken/minikraken_20171019_4GB.tgz"
fi

# All other databases will need to be hosted somehwere before being able to be checked/updated. Currently they are included in the Docker image

# # ANI sketch file / aniDB (150MBs/1.3Gbs)
# if [[ ! -d "${local_DBs}/ANI" ]]; then
# 	cp -r /container_DBs/ANI ${local_DBs}
# fi

# star (6 Mbs)
if [[ ! -d "${local_DBs}/star" ]]; then
	cp -r /container_DBs/star ${local_DBs}
fi

# custom singularity images (3.6GBs)
if [[ ! -d "${local_DBs}/custom_singularities" ]]; then
	cp -r /container_DBs/custom_singularities ${local_DBs}
fi

if [[ ! -d "${local_DBs}/MMB_Bugs.txt" ]]; then
	cp -r /container_DBs/MMB_Bugs.txt ${local_DBs}
fi

if [[ ! -d "${local_DBs}/taxes.csv" ]]; then
	cp -r /container_DBs/taxes.csv ${local_DBs}
fi

if [[ ! -d "${local_DBs}/phiX.fasta" ]]; then
	cp -r /container_DBs/phiX.fasta ${local_DBs}
fi

if [[ ! -d "${local_DBs}/adapters.fasta" ]]; then
	cp -r /container_DBs/adapters.fasta ${local_DBs}
fi

# individual Files (21 KBs)
		#phiX (6 KBs)
		#Adapters (6KBs)
		#taxes (3KBs)
		#bugs (6KBs)
