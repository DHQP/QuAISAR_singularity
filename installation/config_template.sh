#! /usr/bin/env bash
#!/bin/bash -l

#
# Description: Script to consolidate all configuration type settings for QuAISAR pipeline and any tools contained within
# 	Just needs to be sourced within a script to acquire all variables stored within
#
# Usage: . ./config.sh
#
# Output location: No output created
#
# Modules required:
#
# v1.0.4 (3/30/2020)
#
# Created by Nick Vlachos (nvx4@cdc.gov)
#

# Get hostname to help determine if certain tools can be run and how to specifically get others to run with the right options
hostname=$(hostname -f)
host=$(echo "${hostname}" | cut -d'.' -f1)


############# General Options #############
#shortcut to processed samples folder
output_dir=""
# Locations of all scripts and necessary accessory files
src="$(pwd)"
# Local databases that are necessary for pipeline...ANI, BUSCO, star, adapters, phiX
local_DBs=""
# Number of processors requested by numerous applications within the pipeline
procs=4 # Number of processors



############# Application Specific Options #############


#####BBDUK specific config options #####
#requested memory size block
bbduk_mem=Xmx2g
#Kmer length (1-31). Larger Kmer size results in greater specificity
bbduk_k=31
#hamming distance
bbduk_hdist=1


#####Trimmomatic specific config options #####
#Which scoring scale to use
trim_phred="phred33"
#Location of the Adapter FASTA file used for trimming
trim_adapter_location="${local_DBs}/adapters.fasta"
#Seeding mismatches
trim_seed_mismatch=2
#palindrome clip threshold
trim_palindrome_clip_threshold=30
#simple clip threshold
trim_simple_clip_threshold=10
#Minimum adapter length (in palindrome)
trim_min_adapt_length=8
#Keeps read if complete forward and reverse sequences are palindromes
trim_complete_palindrome=TRUE
#Window Size
trim_window_size=20
#Window quality
trim_window_qual=30
#Specifies minimum quality to keep a read
trim_leading=20
#Specifies minimum quality to keep a read
trim_trailing=20
#Specifies minimum length to keep a read
trim_min_length=50


##### SPAdes specific options #####
#Phred quality offset based on scoring used
spades_phred_offset=33
#Coverage threshold (positive float, off or auto)
spades_cov_cutoff="auto"
#Max memory in Gbs
spades_max_memory=8

##### ANI specific options #####
#Max number of samples to be kept (not including source sample) when creating the mash tree
max_ani_samples=20
ani_coverage_threshold=70
# Shortcuts used to find newest REFSEQ mash sketch
REFSEQ=$(find ${local_DBs}/ANI/REFSEQ_*.msh -maxdepth 1 -type f -printf '%p\n' | sort -k2,2 -rt '_' -n | head -n 1)
REFSEQ_date=$(echo "${REFSEQ}" | rev | cut -d'/' -f1 | rev | cut -d'_' -f2 | cut -d'.' -f1)

##### c-SSTAR standard settings #####
# gapped versus ungapped
csstar_gapping="gapped"
# Value used for %id cutoff in csstar (Identity % 100(p), 99(u), 98(h), 95(m), 80(low))
csim=98


##### kraken unclassified threshold ######
# Will throw a warning flag during run summary if percent of unclassified reads/contigs are above this value
unclass_flag=30
# Will throw a warning flag during run summary if percent of 2nd organism is above this value
contamination_threshold=25
# MiniKraken DB (smaller, but effective option)
kraken_DB_path="${local_DBs}/kraken/minikraken_20171019_4GB"
kraken_DB=$(echo "${kraken_DB_path}" | rev | cut -d'/' -f1 | rev)


##### plasmidFinder ######
#percent identity to match
plasmidFinder_identity=.95
#percent length minimum (although not found in command line version, yet)
plasmidFinder_length=60


# Shortcuts used to reference NEWEST AR database
ResGANNCBI_srst2=$(find ${local_DBs}/star/ResGANNCBI_*_srst2.fasta -maxdepth 1 -type f -printf '%p\n' | sort -k2,2 -rt '_' -n | head -n 1)
ResGANNCBI_srst2_filename=$(echo "${ResGANNCBI_srst2}" | rev | cut -d'/' -f1 | rev | cut -d'_' -f1,2)
