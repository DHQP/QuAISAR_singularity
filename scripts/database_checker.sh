#!/bin/bash -l

#$ -o database_checker.out
#$ -e database_checker.err
#$ -N database_checker
#$ -cwd
#$ -q short.q

#
# Description: Script checks for all databases used by QuAISAR pipeline and sets up any missing ones
#
# Usage ./database_checker.sh path_to_database_folder [-i]
#	optional -i is for chencking AND installing, otherwise script just checks for existence..no downloading
#
# Modules required: None
#
# v1.0.1 (05/08/2020)
#
# Created by Nick Vlachos (nvx4@cdc.gov)
#

do_download="false"
# Checks for proper argumentation
if [[ $# -eq 0 ]]; then
	echo "No argument supplied to $0, exiting"
	exit 113
elif [[ "${1}" = "-h" ]]; then
	echo "Usage ./database_checker.sh path_to_config_file [-i]"
	echo "-i is too install databases, otherwise script just checks for existence"
	exit 0
elif [[ ! -d "${1}" ]]; then
	echo "No folder exists as ${1}...exiting"
	exit 113
else
	path_to_DBs="${1}"
fi

if [[ "${2}" == "-i" ]]; then
	do_download="true"
fi

# Shows where databases should be (installed)
echo "${path_to_DBs}"
missing_DBS=()

# Check for parent directory
if [[ ! -d ${path_to_DBs} ]]; then
	mkdir -p ${path_to_DBs}
fi

# # # # Check for BUSCO
#busco_taxa=(bacteria_odb10.2019-06-26 alteromonadales_odb10.2019-04-24 bacillales_odb10.2019-04-24 bacilli_odb10.2019-04-24 bacteroidetes_odb10.2019-04-24 betaproteobacteria_odb10.2019-04-24 burkholderiales_odb10.2019-04-24 campylobacterales_odb10.2019-04-24 clostridiales_odb10.2019-04-24 clostridia_odb10.2019-04-24 corynebacteriales_odb10.2019-04-24 enterobacterales_odb10.2019-04-24 epsilonproteobacteria_odb10.2019-04-24 firmicutes_odb10.2019-04-24 flavobacteriales_odb10.2019-04-24  #flavobacteriia_odb10.2019-04-24 gammaproteobacteria_odb10.2019-04-24 lactobacillales_odb10.2019-04-24 neisseriales_odb10.2019-04-24 proteobacteria_odb10.2019-04-24 pseudomonadales_odb10.2019-04-24 xanthomonadales_odb10.2019-04-24 actinobacteria_class_odb10.2019-04-24 )

# Uploaded 10-16-2020
busco_updated_on1="2020-03-06"
busco_taxa=(acidobacteria_odb10.${busco_updated_on1} actinobacteria_phylum_odb10.${busco_updated_on1} actinobacteria_class_odb10.${busco_updated_on1} bacteroidetes_odb10.${busco_updated_on1} chlamydiae_odb10.${busco_updated_on1} chlorobi_odb10.${busco_updated_on1} chloroflexi_odb10.${busco_updated_on1} cyanobacteria_odb10.${busco_updated_on1} firmicutes_odb10.${busco_updated_on1} fusobacteria_odb10.${busco_updated_on1} planctomycetes_odb10.${busco_updated_on1} proteobacteria_odb10.${busco_updated_on1} spirochaetes_odb10.${busco_updated_on1} synergistetes_odb10.${busco_updated_on1} tenericutes_odb10.${busco_updated_on1} thermotogae_odb10.${busco_updated_on1} verrucomicrobia_odb10.${busco_updated_on1} alphaproteobacteria_odb10.${busco_updated_on1} aquificae_odb10.${busco_updated_on1} bacilli_odb10.${busco_updated_on1} bacteroidia_odb10.${busco_updated_on1} betaproteobacteria_odb10.${busco_updated_on1} clostridia_odb10.${busco_updated_on1} coriobacteriia_odb10.${busco_updated_on1} cytophagia_odb10.${busco_updated_on1} deltaproteobacteria_odb10.${busco_updated_on1} epsilonproteobacteria_odb10.${busco_updated_on1} flavobacteriia_odb10.${busco_updated_on1} gammaproteobacteria_odb10.${busco_updated_on1} mollicutes_odb10.${busco_updated_on1} sphingobacteriia_odb10.${busco_updated_on1} spirochaetia_odb10.${busco_updated_on1} tissierellia_odb10.${busco_updated_on1} alteromonadales_odb10.${busco_updated_on1} bacillales_odb10.${busco_updated_on1} bacteroidales_odb10.${busco_updated_on1} burkholderiales_odb10.${busco_updated_on1} campylobacterales_odb10.${busco_updated_on1} cellvibrionales_odb10.${busco_updated_on1} chromatiales_odb10.${busco_updated_on1} chroococcales_odb10.${busco_updated_on1} clostridiales_odb10.${busco_updated_on1} coriobacteriales_odb10.${busco_updated_on1} corynebacteriales_odb10.${busco_updated_on1} cytophagales_odb10.${busco_updated_on1} desulfobacterales_odb10.${busco_updated_on1} desulfovibrionales_odb10.${busco_updated_on1} desulfuromonadales_odb10.${busco_updated_on1} enterobacterales_odb10.${busco_updated_on1} entomoplasmatales_odb10.${busco_updated_on1} flavobacteriales_odb10.${busco_updated_on1} fusobacteriales_odb10.${busco_updated_on1} lactobacillales_odb10.${busco_updated_on1} legionellales_odb10.${busco_updated_on1} micrococcales_odb10.${busco_updated_on1} mycoplasmatales_odb10.${busco_updated_on1} neisseriales_odb10.${busco_updated_on1} nitrosomonadales_odb10.${busco_updated_on1} nostocales_odb10.${busco_updated_on1} oceanospirillales_odb10.${busco_updated_on1} oscillatoriales_odb10.${busco_updated_on1} pasteurellales_odb10.${busco_updated_on1} propionibacteriales_odb10.${busco_updated_on1} pseudomonadales_odb10.${busco_updated_on1} rhizobiales_odb10.${busco_updated_on1} rhodobacterales_odb10.${busco_updated_on1} rhodospirillales_odb10.${busco_updated_on1} rickettsiales_odb10.${busco_updated_on1} selenomonadales_odb10.${busco_updated_on1} sphingomonadales_odb10.${busco_updated_on1} )



#echo "${#busco_taxa[@]}"

for odb_info in "${busco_taxa[@]}"; do
	# Check for top level bacteria database
	#echo ${odb_info}
	taxa=$(echo "$odb_info" | cut -d'.' -f1)
	db_date=$(echo "$odb_info" | cut -d'.' -f2)
	if [[ ! -d "${path_to_DBs}/BUSCO/${taxa}" ]]; then
		if [[ "${do_download}" = "true" ]]; then
			if [[ ! -d "${path_to_DBs}/BUSCO" ]]; then
				mkdir "${path_to_DBs}/BUSCO"
			fi
			cd "${path_to_DBs}/BUSCO"
			echo "Downloading latest BUSCO database for ${taxa} (wget http://busco-data.ezlab.org/v4/data/lineages/${taxa}.${db_date}.tar.gz)"
			wget "http://busco-data.ezlab.org/v4/data/lineages/${taxa}.${db_date}.tar.gz"
			# Dont know how to handle this one outlier (only one to specify a level in the filename) - ALl OUR bugs are in class Actinobacteria too
		else
			echo "Missing latest BUSCO database for ${taxa}"
			missing_DBS=("${missing_DBS[@]}" "BUSCO-${taxa}")
		fi
	else
		echo "BUSCO has latest ${taxa}_odb10 as of 3/15/2020"
	fi
done

# Unpack and cleanup
if [[ "${do_download}" = "true" ]]; then
	find ${path_to_DBs}/BUSCO/ -name '*.gz' -exec tar xzf {} \;
	#mv ${path_to_DBs}/BUSCO/actinobacteria_class_odb10 ${path_to_DBs}/BUSCO/actinobacteria_odb10
	find ${path_to_DBs}/BUSCO/ -name '*.gz' -exec rm {} \;
fi

# All other databases will need to be hosted somehwere before being able to be checked/updated. Currently they are included in the Docker image

# Test index 0 is filename output, 1 is MEGA link, 2 is onedrive link, 3 is drive link
link_index=4
# Lists of links to test for downloading (MEGA OneDrive Google)
###bbtools_links=("bbtools.simg" "https://mega.nz/file/0r5UCYIR#zn3LHj7RHKAMR-VkDGSc-5lUmWaE12A3jBPQOCJaZOk" "https://onedrive.live.com/download?cid=89BB0F0D841B2A3B&resid=89BB0F0D841B2A3B%21106&authkey=AHXpg4F2NHk28Vw" "https://drive.google.com/uc?export=download&id=1QCvz1LRidSmeXhzMrfm7GtEDff4ik3zs" "ftp://ftp.cdc.gov/pub/QUAISAR-FTP/bbtools.simg")
###blast_links=("blast-2.9.0-docker.img" "https://mega.nz/file/gyhCVIQR#1n-m6DEI1LA6HOiEE40i9x3fv5iXYFZsWT9sKfsNs_M" "https://onedrive.live.com/download?cid=89BB0F0D841B2A3B&resid=89BB0F0D841B2A3B%21111&authkey=ADVSQ-oAmV3VAJk" "https://drive.google.com/uc?export=download&id=1-Ic9CxvcR4ubNsJX6k2XubWKhFXS_GVc" "ftp://ftp.cdc.gov/pub/QUAISAR-FTP/blast-2.9.0-docker.img")
###bowtie2_links=("bowtie2-2.2.9-biocontainers.simg" "https://mega.nz/file/92hgXAQJ#XThfqohBWcpD3kRzgUv4RjscDnmS2Xl5lNtBHVnvNuw" "https://onedrive.live.com/download?cid=89BB0F0D841B2A3B&resid=89BB0F0D841B2A3B%21113&authkey=AKpz8agXXcDjON0" "https://drive.google.com/uc?export=download&id=1rQncOI9oqQKRDFLvxJFrZSZh2We5WA_S" "ftp://ftp.cdc.gov/pub/QUAISAR-FTP/bowtie2-2.2.9-biocontainers.simg")
cSSTAR_links=("csstar.simg" "https://mega.nz/file/12wkUSDB#huoDBxj6keneY9h0hehwBWPoZ_n5zTpOiELePL0szFs" "https://onedrive.live.com/download?cid=89BB0F0D841B2A3B&resid=89BB0F0D841B2A3B%21117&authkey=AFzS-KUAFnQdSnU" "https://drive.google.com/uc?export=download&id=1b--ivEqsFdnPPB85Ds340J8Pu3BzSRoI" "ftp://ftp.cdc.gov/pub/QUAISAR-FTP/csstar.simg")
entrez_links=("entrez_taxon.simg" "https://mega.nz/file/gy4RQQrY#JZWvV4-PbeMfOJjnRj6qjZ9jzkXZGFgLbURGPyQu42E" "https://onedrive.live.com/download?cid=89BB0F0D841B2A3B&resid=89BB0F0D841B2A3B%21116&authkey=AGPABuxJPVYnJgA" "https://drive.google.com/uc?export=download&id=1Earz_6jrnTkIjV_oLTmnuY9P2NSMY614" "ftp://ftp.cdc.gov/pub/QUAISAR-FTP/entrez_taxon.simg")
###GAMA_links=("GAMA_quaisar.simg" "https://mega.nz/file/R3wmFQJb#yY3gQ1tFvIPxeKSEUydezyTh5fnVANBOA0LV7dmHHFk" "https://onedrive.live.com/download?cid=89BB0F0D841B2A3B&resid=89BB0F0D841B2A3B%21107&authkey=AJKsAgfit5oEhCs" "https://drive.google.com/uc?export=download&id=1kzEloxzGIVPaV-1PngLt9FAn2dpGeguj" "ftp://ftp.cdc.gov/pub/QUAISAR-FTP/GAMA_quaisar.simg")
gottcha_links=("gottcha.simg" "https://mega.nz/file/EyxEDCrC#Q2kkGwzDB0HdLL3q9U2uRf7gd1orHbFK_voCCTIBErc" "https://onedrive.live.com/download?cid=89BB0F0D841B2A3B&resid=89BB0F0D841B2A3B%21108&authkey=AMIQfW5e3Yi-sH4" "https://drive.google.com/uc?export=download&id=1SnGgM60JUdhe0y5EoZprcPpG326pG79j" "ftp://ftp.cdc.gov/pub/QUAISAR-FTP/gottcha.simg")
plasmidFinder_links=("plasmidFinder_with_DB.simg" "https://mega.nz/file/kugiyYZK#um_iss6jLcs4P3_qL7M5EYHICcyYJz0cHyCmUaR4ovg" "https://onedrive.live.com/download?cid=89BB0F0D841B2A3B&resid=89BB0F0D841B2A3B%21112&authkey=AHORh2N541BayFw" "https://drive.google.com/uc?export=download&id=1wZAzkiTD2rLxkgDon16rZErMjg986yit" "ftp://ftp.cdc.gov/pub/QUAISAR-FTP/plasmidFinder_with_DB.simg")
###QUAST_links=("QUAST5.simg" "https://mega.nz/file/8rw0kQxL#1p-zUtABJb9sLmwkeAojSMmFJ8oRkZaOtVinT0Jo1NY" "https://onedrive.live.com/download?cid=89BB0F0D841B2A3B&resid=89BB0F0D841B2A3B%21110&authkey=AM_c7ih_fP4JxrE" "https://drive.google.com/uc?export=download&id=1JhxjA2xt4dsjpO96egc7Iyz9wOEbnXXF" "ftp://ftp.cdc.gov/pub/QUAISAR-FTP/QUAST5.simg")
###srst2_links=("srst2.simg" "https://mega.nz/file/Y6hg3CCb#6lLqih6Dv5AYOs0hfJiBZD7BkxR8k4wwhTEkJKKmwls" "https://onedrive.live.com/download?cid=89BB0F0D841B2A3B&resid=89BB0F0D841B2A3B%21109&authkey=AINwP6LEwO1bDgI" "https://drive.google.com/uc?export=download&id=1Kobw285kXNy7yxHHxhRxYY2PjjEmVtYK" "ftp://ftp.cdc.gov/pub/QUAISAR-FTP/srst2.simg")
srst2_links=("srst2.simg" " " " " " " "ftp://ftp.cdc.gov/pub/QUAISAR-FTP/srst2_corrected.sif")
ANI_links=("REFSEQ_ANI.msh.gz" "https://mega.nz/file/puxixKaT#tUbaDQ1YV2TpxgpHyhlOI1ryTfaBP7RhBgD9_Psimhc" "https://onedrive.live.com/download?cid=89BB0F0D841B2A3B&resid=89BB0F0D841B2A3B%21114&authkey=AB5T8jQOePfzxSg" "https://drive.google.com/uc?export=download&id=161jVEG-AV38qNxIcNHSM-T0hn1RoyIv0" "ftp://ftp.cdc.gov/pub/QUAISAR-FTP/REFSEQ_ANI.msh.gz")
pubmlst_links=("pubmlst.tar.gz" "https://mega.nz/file/M2h2mYZT#SJ4ohNn60WsdHovWxNKp6sQTtcA5tfk6WY-iECA9zEw" "https://onedrive.live.com/download?cid=89BB0F0D841B2A3B&resid=89BB0F0D841B2A3B%21115&authkey=AGCIPp4ZdSRdGHc" "https://drive.google.com/uc?export=download&id=1DoqUliXXJSWEsZFCoGSnakzjySbywgx2" "ftp://ftp.cdc.gov/pub/QUAISAR-FTP/pubmlst.tar.gz")
sstar_links=("sstar.tar.gz" "https://mega.nz/file/Z7gEWArZ#MfOJld0JsjtYMXI7vzkr-N2f8oKCpTgM1zYobL6fX3E" "https://onedrive.live.com/download?cid=89BB0F0D841B2A3B&resid=89BB0F0D841B2A3B%21122&authkey=AG1SomRvYC1gMxM" "https://drive.google.com/uc?export=download&id=1WXqL4bdT-eO_zyIk-csLPdZatj-JTyPM" "ftp://ftp.cdc.gov/pub/QUAISAR-FTP/sstar.tar.gz")
adapters_links=("adapters.fasta" "https://mega.nz/file/AqwyXA7S#Ao8VR1JELCeos6ISbDE3e1r7LNXSplfA5pc5m8ICb5w" "https://onedrive.live.com/download?cid=89BB0F0D841B2A3B&resid=89BB0F0D841B2A3B%21121&authkey=AO1Xn0MaUceaHVw" "https://drive.google.com/uc?export=download&id=1Ec_tQoL-fsMYJArP92-1g--SHoEiEhIY" "ftp://ftp.cdc.gov/pub/QUAISAR-FTP/adapters.fasta")
#MMBbugs_links=("MMB_Bugs.txt" "https://mega.nz/file/ArxWwCbD#gtCCSQAKZCKB_euY2tryANlq_R4hZkN_HdZhEfPjo1k" "https://onedrive.live.com/download?cid=89BB0F0D841B2A3B&resid=89BB0F0D841B2A3B%21128&authkey=ADtwVO9h0vMIU4M" "https://drive.google.com/uc?export=download&id=1RYL8timv8CLhvIyXgRKiTUWxwerrvEcf" "ftp://ftp.cdc.gov/pub/QUAISAR-FTP/MMB_Bugs.txt")
phiX_links=("phiX.fasta" "https://mega.nz/file/drg00AIK#rH95tA3qTuE7SDGvcUIDuGfENm4un4MPYyI5G1tmCB4" "https://onedrive.live.com/download?cid=89BB0F0D841B2A3B&resid=89BB0F0D841B2A3B%21119&authkey=AEYBttLLp5mkKM0" "https://drive.google.com/uc?export=download&id=111clfbuw7sjUoQN0Uw8NVhZINJAf739G" "ftp://ftp.cdc.gov/pub/QUAISAR-FTP/phiX.fasta")
taxes_links=("taxes.csv" "https://mega.nz/file/1v5UFIgY#K5pFtVNLrKP5kB6BH5S_f_nqwVNvLRhnW7xkd8dSteo" "https://onedrive.live.com/download?cid=89BB0F0D841B2A3B&resid=89BB0F0D841B2A3B%21127&authkey=AOl2_e_i2ubfrt4" "https://drive.google.com/uc?export=download&id=1-lnexzeEquhcvJT-UE3wx8YD804rr5ts" "ftp://ftp.cdc.gov/pub/QUAISAR-FTP/taxes.csv")
ratio_links=("NCBI_Assembly_stats.txt.gz" "not_saved_to_mega" "not_saved_to_onedrive" "not_saved_to_google" "ftp://ftp.cdc.gov/pub/QUAISAR-FTP/NCBI_Assembly_stats.txt.gz")

wget_options=""
if [[ ${link_index} -eq 1 ]]; then
	echo "Unknown how to process"
elif [[ ${link_index} -ge 2 ]]; then
	wget_options="--no-check-certificate"
elif [[ ${link_index} -eq 3 ]]; then
	echo "None needed for google drive"
elif [[ ${link_index} -eq 4 ]]; then
	echo "None needed for CDC FTP"
fi

# star (6 Mbs)
if [[ ! -d "${path_to_DBs}/star" ]]; then
	#cp -r /container_DBs/star ${path_to_DBs}
	if [[ "${do_download}" = "true" ]]; then
		echo "Copying latest NAR-AR database"
		cd ${path_to_DBs}
		wget "${wget_options}" -O "${sstar_links[0]}" "${sstar_links[${link_index}]}"
		tar -zxvf sstar.tar.gz
		rm sstar.tar.gz
	else
		echo "Missing latest NAR-AR database"
		missing_DBS=("${missing_DBS[@]}" "NAR-AR")
	fi
else
	echo "NAR-AR database installed"
fi

# if [[ ! -f "${path_to_DBs}/MMB_Bugs.txt" ]]; then
# 	#cp -r /container_DBs/MMB_Bugs.txt ${path_to_DBs}
# 	if [[ "${do_download}" = "true" ]]; then
# 		echo "Copying MMB_bugs"
# 		cd ${path_to_DBs}
# 		wget "${wget_options}" -O "${MMBbugs_links[0]}" "${MMBbugs_links[${link_index}]}"
# 	else
# 		echo "Missing MMB_Bugs"
# 		missing_DBS=("${missing_DBS[@]}" "MMB_Bugs")
# 	fi
# else
# 	echo "MMB_Bugs installed"
# fi

if [[ ! -f "${path_to_DBs}/taxes.csv" ]]; then
	#cp -r /container_DBs/taxes.csv ${path_to_DBs}
	if [[ "${do_download}" = "true" ]]; then
		echo "Copying taxes"
		cd ${path_to_DBs}
		wget "${wget_options}" -O "${taxes_links[0]}" "${taxes_links[${link_index}]}"
	else
		echo "Missing taxes"
		missing_DBS=("${missing_DBS[@]}" "taxes")
	fi
else
	echo "taxes installed"
fi

if [[ ! -f "${path_to_DBs}/phiX.fasta" ]]; then
	#cp -r /container_DBs/phiX.fasta ${path_to_DBs}
	if [[ "${do_download}" = "true" ]]; then
		echo "Copying phiX.fasta"
		cd ${path_to_DBs}
		wget "${wget_options}" -O "${phiX_links[0]}" "${phiX_links[${link_index}]}"
	else
		echo "Missing phiX"
		missing_DBS=("${missing_DBS[@]}" "phiX")
	fi
else
	echo "phiX installed"
fi

if [[ ! -f "${path_to_DBs}/adapters.fasta" ]]; then
	#cp -r /container_DBs/adapters.fasta ${path_to_DBs}
	if [[ "${do_download}" = "true" ]]; then
		echo "Copying adapters.fasta"
		cd ${path_to_DBs}
		wget "${wget_options}" -O "${adapters_links[0]}" "${adapters_links[${link_index}]}"
	else
		echo "Missing adapters"
		missing_DBS=("${missing_DBS[@]}" "adapters")
	fi
else
	echo "adapters installed"
fi

if [[ ! -d "${path_to_DBs}/ratio_DBs" ]]; then
	#cp -r /container_DBs/adapters.fasta ${path_to_DBs}
	if [[ "${do_download}" = "true" ]]; then
		echo "Copying NCBI_Assembly_stats"
		mkdir -p ${path_to_DBs}/ratio_DBs
		cd ${path_to_DBs}/ratio_DBs
		wget "${wget_options}" -O "${ratio_links[0]}" "${ratio_links[${link_index}]}"
		gunzip -N NCBI_Assembly_stats.txt.gz
	else
		echo "Missing ratios file"
		missing_DBS=("${missing_DBS[@]}" "NCBI_ratio")
	fi
else
	echo "ratios installed"
fi

if [[ ! -d "${path_to_DBs}/ANI" ]]; then
	#cp -r /container_DBs/ANI ${path_to_DBs}
	if [[ "${do_download}" = "true" ]]; then
		echo "Copying latest REFSEQ sketch database (ANI)"
		mkdir ${path_to_DBs}/ANI
		cd ${path_to_DBs}/ANI
		if [[ ${link_index} -eq 3 ]]; then
			query=`curl -k -c ./cookie.txt -s -L "${ANI_links[3]}" \
			| perl -nE'say/uc-download-link.*? href="(.*?)\">/' \
			| sed -e 's/amp;//g' | sed -n 2p`
			url="https://drive.google.com$query"
			curl -k -b ./cookie.txt -L -o ${ANI_links[0]} $url
			rm ./cookie.txt
		else
			wget "${wget_options}" -O "${ANI_links[0]}" "${ANI_links[${link_index}]}"
		fi
		gunzip -N REFSEQ_ANI.msh.gz
	else
		echo "Missing latest REFSEQ sketch database (ANI)"
		missing_DBS=("${missing_DBS[@]}" "REFSEQ-ANI")
	fi
else
	echo "ANI REFSEQ sketch database installed"
fi

if [[ ! -d "${path_to_DBs}/pubmlst" ]]; then
	#cp -r /container_DBs/pubmlst ${path_to_DBs}
	if [[ "${do_download}" = "true" ]]; then
		echo "Copying pubMLST"
		cd ${path_to_DBs}
		wget {wget_options} -O "${pubmlst_links[0]}" "${pubmlst_links[${link_index}]}"
		tar -zxvf pubmlst.tar.gz
		rm pubmlst.tar.gz
		# for scheme in pubmlst/*/*.gz; do
		# 	gunzip $scheme
		# done
	else
		echo "Missing pubMLST"
		missing_DBS=("${missing_DBS[@]}" "pubMLST")
	fi
else
	echo "pubMLST installed"
fi

#singularities=(bbtools.simg+${bbtools_links[${link_index}]}+76 blast-2.9.0-docker.img+${blast_links[${link_index}]}+94 bowtie2-2.2.9-biocontainers.simg+${bowtie2_links[${link_index}]}+364 cSSTAR.simg+${cSSTAR_links[${link_index}]}+688 entrez_taxon.simg+${entrez_links[${link_index}]}+239 GAMA_quaisar.simg+${GAMA_links[${link_index}]}+242 gottcha.simg+${gottcha_links[${link_index}]}+208 plasmidFinder_with_DB.simg+${plasmidFinder_links[${link_index}]}+805 QUAST5.simg+${QUAST_links[${link_index}]}+345 srst2.simg+${srst2_links[${link_index}]}+262)

singularities=(cSSTAR.simg+${cSSTAR_links[${link_index}]}+688 entrez_taxon.simg+${entrez_links[${link_index}]}+239 plasmidFinder_with_DB.simg+${plasmidFinder_links[${link_index}]}+805 srst2.simg+${srst2_links[${link_index}]}+262) #gottcha.simg+${gottcha_links[${link_index}]}+208

for simage_info in "${singularities[@]}"; do
	# custom singularity images (3.6GBs)
	simage=$(echo "${simage_info}" | cut -d'+' -f1)
	url_link=$(echo "${simage_info}" | cut -d'+' -f2)
	size=$(echo "${simage_info}" | cut -d'+' -f3)
	echo -e "${simage}\n${url_link}\n${size}\n"

	if [[ ! -f "${path_to_DBs}/singularities/${simage}" ]]; then
		#cp -r /container_DBs/custom_singularities ${path_to_DBs}
		if [[ "${do_download}" = "true" ]]; then
			if [[ ! -d "${path_to_DBs}/singularities" ]]; then
				mkdir "${path_to_DBs}/singularities"
				cd "${path_to_DBs}/singularities"
			fi
			echo "Copying custom singularity image ${simage}"
			if [[ ${size} -ge 100 ]] && [[ ${link_index} -eq 3 ]]; then
				echo "Too big, special command"
				#wget --save-cookies cookies.txt "${url_link}" -O- | sed -rn 's/.*confirm=([0-9A-Za-z_]+).*/\1/p' > confirm.txt
		 		#wget --load-cookies cookies.txt -O ${simage} '${url_link}'&'confirm='$(<confirm.txt)

				query=`curl -k -c ./cookie.txt -s -L "${url_link}" \
				| perl -nE'say/uc-download-link.*? href="(.*?)\">/' \
				| sed -e 's/amp;//g' | sed -n 2p`
				url="https://drive.google.com$query"
				curl -k -b ./cookie.txt -L -o ${simage} $url
				rm ./cookie.txt
			else
				echo "Normal command -just testing"
				wget "${wget_options}" -O ${simage} "${url_link}"
			fi
		else
			echo "Missing custom singularity image ${simage}"
			missing_DBS=("${missing_DBS[@]}" "singularities-${simage}")
		fi
	else
		echo "custom singularity image ${simage} installed"
	fi
done

# Check to see if kraken mini database is installed
if [[ ! -d "${path_to_DBs}/kraken" ]]; then
	if [[ "${do_download}" = "true" ]]; then
		mkdir "${path_to_DBs}/kraken"
		cd "${path_to_DBs}/kraken"
		echo "Downloading latest (mini)kraken database (wget ftp://ftp.cdc.gov/pub/QUAISAR-FTP/minikraken_20171019_4GB.tgz)"

		wget "ftp://ftp.cdc.gov/pub/QUAISAR-FTP/minikraken_20171019_4GB.tgz"
		if [[ ! -f "minikraken_20171019_4GB.tgz" ]]; then
			curl -k -O "ftp://ftp.cdc.gov/pub/QUAISAR-FTP/minikraken_20171019_4GB.tgz"
		fi
		tar xzf minikraken_20171019_4GB.tgz

		mv minikraken_20171013_4GB minikraken_20171019_4GB

		rm minikraken_20171019_4GB.tgz
	else
		echo "Missing latest kraken database"
		missing_DBS=("${missing_DBS[@]}" "kraken")
	fi
else
	echo "kraken database is installed"
fi

# ##### Currently down.....and has been a while
# # Check to see if gottcha database is installed
# if [[ ! -d "${path_to_DBs}/gottcha" ]]; then
# 	if [[ "${do_download}" = "true" ]]; then
# 		cd "${path_to_DBs}"
# 		# Original LANL hosted address that has been down a good while
# 	 	#wget -P "${path_to_DBs}/gottcha" "https://edge-dl.lanl.gov/gottcha/GOTTCHA_database_v20150825/GOTTCHA_BACTERIA_c4937_k24_u30_xHUMAN3x.species.tar.gz"
# 		# Temporary mirror until original is fixed
# 		echo "Downloading latest gottcha database (wget ftp://ftp.cdc.gov/pub/QUAISAR-FTP/gottcha_bac_arc_v1.tar.gz)"
# 		wget "ftp://ftp.cdc.gov/pub/QUAISAR-FTP/gottcha_bac_arc_v1.tar.gz"
# 		tar xzf gottcha_bac_arc_v1.tar.gz
# 		rm gottcha_bac_arc_v1.tar.gz
# 		#mv gottcha/gottcha_db ./
# 		#rm -r gottcha
# 		#mv gottcha_db gottcha
# 		#rm gottcha.dbprofile.out
# 		# Need to find sa place to host genus_Lookup.tar.gz
# 	else
# 		echo "Missing gottcha database"
# 		missing_DBS=("${missing_DBS[@]}" "gottcha")
# 	fi
# else
# 	echo "gottcha database installed"
# fi

ls ${path_to_DBs}
#chmod -Rx ${path_to_DBs}/*

echo "There are ${#missing_DBS[@]} missing databases (${missing_DBS[@]})"
