#!/bin/bash -l

 #
 # Description: Script to clean extraneous info from a SPAdes run
 #
 # Usage: ./gulvic_SPAdes_cleaner.sh directory_to_clean (typically Project/Sample_name/Assembly within our structure)
 #
 # Output location: No output created
 #
 # Modules required: None
 #
 # v1.0 (10/3/2019)
 #
 # Created by Chris Gulvik
 #

if [[ "$1" == "" || "$1" == "--help" || "$1" == "-h" ]]; then
 	echo "
 	Usage: `basename $0` input_dir

 	Given a path containing SPAdes assemblies, recursively
 	searches within the input dir to remove most dirs and files
 	SPAdes creates but maintains essential log files, contigs,
 	and scaffolds FastA files."
	exit 1
fi

mapfile -t asm_dirs < <(find -L "$1" -type d -regextype posix-extended -regex '\/.*\/K[1-9][0-9]{1,2}' -print | xargs dirname | uniq)

for directory in "${asm_dirs[@]}"; do
	rm -v "$directory"/{before_rr.fasta,contigs.paths,dataset.info,input_dataset.yaml,scaffolds.paths};
done

find -L "$1" -type d -regextype posix-extended -regex '\/.*\/(K[1-9][0-9]{1,2}|misc|tmp)' -print | xargs rm -rv
