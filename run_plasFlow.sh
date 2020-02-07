#!/bin/sh -l

#$ -o pFlow.out
#$ -e pFlow.err
#$ -N pFlow
#$ -cwd
#$ -q short.q

# Import the config file with shortcuts and settings
if [[ ! -f "./config.sh" ]]; then
	cp ./config_template.sh ./config.sh
fi
. ./config.sh


#
# Description: Will attempt to find any plasmids in sample using plasFlow methods
#
# Usage: ./run_plasFlow.sh sample_name run_ID
#
# Output location: default_config.sh_output_location/run_ID/sample_name/plasFlow/
#
# Modules required: PlasFlow/1.1
#
# v1.0.1 (10/9/2019)
#
# Created by Nick Vlachos (nvx4@cdc.gov)
#

# Checks for proper argumentation
if [[ $# -eq 0 ]]; then
	echo "No argument supplied to $0, exiting"
	exit 1
elif [[ -z "${1}" ]]; then
	echo "Empty sample name supplied to run_plasmidFlow.sh, exiting"
	exit 1
elif [[ "${1}" = "-h" ]]; then
	echo "Usage is ./run_plasmFlow.sh  sample_name run_ID"
	echo "Output by default is ${output_dir}/miseq_run_ID/sample_name/plasmFlow"
	exit 0
elif [[ -z "${2}" ]]; then
	echo "Empty run_ID supplied to run_plasFlow.sh, exiting"
	exit 1
fi

ml PlasFlow/1.1

# Create output directory
if [[ ! -d "${output_dir}/${2}/${1}/plasFlow" ]]; then
	mkdir "${output_dir}/${2}/${1}/plasFlow"
fi

# Prep any samples that don't have the paired.fq reads
if [ -f "${output_dir}/${2}/${1}/trimmed/${1}_R1_001.paired.fq" ] && [ -f "${output_dir}/${2}/${1}/trimmed/${1}_R2_001.paired.fq" ]; then
	echo "1"
	#gzip -c "${output_dir}/${2}/${1}/trimmed/${1}_R1_001.paired.fq" > "${output_dir}/${2}/${1}/trimmed/${1}_S1_L001_R1_001.fastq.gz"
	#gzip -c "${output_dir}/${2}/${1}/trimmed/${1}_R2_001.paired.fq" > "${output_dir}/${2}/${1}/trimmed/${1}_S1_L001_R2_001.fastq.gz"
elif [ -f "${output_dir}/${2}/${1}/trimmed/${1}_R1_001.paired.fq.gz" ]; then
	echo "2"
	gunzip < "${output_dir}/${2}/${1}/trimmed/${1}_R1_001.paired.fq.gz" > "${output_dir}/${2}/${1}/trimmed/${1}_R1_001.paired.fq"
	if [[ -f "${output_dir}/${2}/${1}/trimmed/${1}_R2_001.paired.fq.gz" ]]; then
	echo "2A"
		gunzip < "${output_dir}/${2}/${1}/trimmed/${1}_R2_001.paired.fq.gz" > "${output_dir}/${2}/${1}/trimmed/${1}_R2_001.paired.fq"
	fi
elif [ -f "${output_dir}/${2}/${1}/trimmed/${1}_R2_001.paired.fq.gz" ]; then
	echo "3"
	gunzip < "${output_dir}/${2}/${1}/trimmed/${1}_R2_001.paired.fq.gz" > "${output_dir}/${2}/${1}/trimmed/${1}_R2_001.paired.fq"
else
	echo "4"
	if [[ -f "${output_dir}/${2}/${1}/FASTQs/${1}_R1_001.fastq.gz" ]] && [[ ! -f "${output_dir}/${2}/${1}/FASTQs/${1}_R1_001.fastq" ]]; then
		gunzip -c "${output_dir}/${2}/${1}/FASTQs/${1}_R1_001.fastq.gz" > "${output_dir}/${2}/${1}/FASTQs/${1}_R1_001.fastq"
	fi
	if [[ -f "${output_dir}/${2}/${1}/FASTQs/${1}_R2_001.fastq.gz" ]] && [[ ! -f "${output_dir}/${2}/${1}/FASTQs/${1}_R2_001.fastq" ]]; then
		gunzip -c "${output_dir}/${2}/${1}/FASTQs/${1}_R2_001.fastq.gz" > "${output_dir}/${2}/${1}/FASTQs/${1}_R2_001.fastq"
	fi
	echo "Running BBDUK and trimmomatic"
	bbduk.sh -"${bbduk_mem}" threads="${procs}" in="${output_dir}/${2}/${1}/FASTQs/${1}_R1_001.fastq" in2="${output_dir}/${2}/${1}/FASTQs/${1}_R2_001.fastq" out="${output_dir}/${2}/${1}/removedAdapters/${1}-noPhiX-R1.fsq" out2="${output_dir}/${2}/${1}/removedAdapters/${1}-noPhiX-R2.fsq" ref="${phiX_location}" k="${bbduk_k}" hdist="${bbduk_hdist}"
	trimmomatic "${trim_endtype}" -"${trim_phred}" -threads "${procs}" "${output_dir}/${2}/${1}/removedAdapters/${1}-noPhiX-R1.fsq" "${output_dir}/${2}/${1}/removedAdapters/${1}-noPhiX-R2.fsq" "${output_dir}/${2}/${1}/trimmed/${1}_R1_001.paired.fq" "${output_dir}/${2}/${1}/trimmed/${1}_R1_001.unpaired.fq" "${output_dir}/${2}/${1}/trimmed/${1}_R2_001.paired.fq" "${output_dir}/${2}/${1}/trimmed/${1}_R2_001.unpaired.fq" ILLUMINACLIP:"${trim_adapter_location}:${trim_seed_mismatch}:${trim_palindrome_clip_threshold}:${trim_simple_clip_threshold}:${trim_min_adapt_length}:${trim_complete_palindrome}" SLIDINGWINDOW:"${trim_window_size}:${trim_window_qual}" LEADING:"${trim_leading}" TRAILING:"${trim_trailing}" MINLEN:"${trim_min_length}"
fi


# Check if sample has original assembly to process plasflow from
if [[ -s "${output_dir}/${2}/${1}/Assembly/${1}_scaffolds_trimmed.fasta" ]]; then
	# Trim contigs a little to 2000 and larger and put through plasflow.
	python3 "${src}/removeShortContigs.py" -i "${output_dir}/${2}/${1}/Assembly/scaffolds.fasta" -t 2000 -s "normal_SPAdes"
	# Renames headers of fasta files accordingly
	python3 "${src}/fasta_headers.py" -i "${output_dir}/${2}/${1}/Assembly/scaffolds.fasta.TRIMMED.fasta" -o "${output_dir}/${2}/${1}/plasFlow/${1}_scaffolds_trimmed_2000.fasta"
	# Removes intermeidate fasta file
	rm -r "${output_dir}/${2}/${1}/Assembly/scaffolds.fasta.TRIMMED.fasta"
	# Run plasflow on newly trimmed assembly file
	PlasFlow.py --input "${output_dir}/${2}/${1}/plasFlow/${1}_scaffolds_trimmed_2000.fasta" --output "${output_dir}/${2}/${1}/plasFlow/${1}_plasFlow_results.tsv" --threshold 0.7
	# Load all necessary modules to complete the realignment portion of analysis

	#module load Python3/3.5.4;
	#module load bowtie2/2.2.9;
	#module load samtools/1.4.1;
	#module load bam2fastq/1.1.0;
	#module load Unicycler/0.4.4;
	#module load gcc/5.5;
	#module load SPAdes/3.13.0;
	#module load racon/1.3.1;
	#module load perl/5.22.1

	ml -Python3/3.5 bowtie2/2.2.9 samtools/1.4.1 bam2fastq/1.1.0 Unicycler/0.4.4 #SPAdes/3.13.0 racon/1.3.1

	mkdir ${output_dir}/${2}/${1}/plasFlow/bowtie2-index/
	bowtie2-build -f "${output_dir}/${2}/${1}/plasFlow/${1}_plasFlow_results.tsv_chromosomes.fasta" "${output_dir}/${2}/${1}/plasFlow/bowtie2-index/bowtie2_${1}_chr"
	mkdir ${output_dir}/${2}/${1}/plasFlow/filtered_reads_70/
	bowtie2 -x "${output_dir}/${2}/${1}/plasFlow/bowtie2-index/bowtie2_${1}_chr" -1 "${output_dir}/${2}/${1}/trimmed/${1}_R1_001.paired.fq" -2 "${output_dir}/${2}/${1}/trimmed/${1}_R2_001.paired.fq" -S "${output_dir}/${2}/${1}/plasFlow/filtered_reads_70/${1}.sam" -p 12 --local
	samtools view -bS "${output_dir}/${2}/${1}/plasFlow/filtered_reads_70/${1}.sam" > "${output_dir}/${2}/${1}/plasFlow/filtered_reads_70/${1}.bam"
	bam2fastq --no-aligned -o "${output_dir}/${2}/${1}/plasFlow/filtered_reads_70/${1}_R#_bacterial.fastq" "${output_dir}/${2}/${1}/plasFlow/filtered_reads_70/${1}.bam"
	unicycler -1 "${output_dir}/${2}/${1}/plasFlow/filtered_reads_70/${1}_R_1_bacterial.fastq" -2 "${output_dir}/${2}/${1}/plasFlow/filtered_reads_70/${1}_R_2_bacterial.fastq" -o "${output_dir}/${2}/${1}/plasFlow/Unicycler_assemblies/${1}_uni_assembly"
	mv "${output_dir}/${2}/${1}/plasFlow/Unicycler_assemblies/${1}_uni_assembly/assembly.fasta" "${output_dir}/${2}/${1}/plasFlow/Unicycler_assemblies/${1}_uni_assembly/${1}_plasmid_assembly_original.fasta"
	mv "${output_dir}/${2}/${1}/plasFlow/Unicycler_assemblies/${1}_uni_assembly/assembly.gfa" "${output_dir}/${2}/${1}/plasFlow/Unicycler_assemblies/${1}_uni_assembly/${1}_assembly.gfa"
	python3 ${src}/fasta_headers_plasFlow.py -i "${output_dir}/${2}/${1}/plasFlow/Unicycler_assemblies/${1}_uni_assembly/${1}_plasmid_assembly_original.fasta" -o "${output_dir}/${2}/${1}/plasFlow/Unicycler_assemblies/${1}_uni_assembly/${1}_plasmid_assembly.fasta"
	python3 "${src}/removeShortContigs.py" -i "${output_dir}/${2}/${1}/plasFlow/Unicycler_assemblies/${1}_uni_assembly/${1}_plasmid_assembly.fasta" -t 500 -s "plasFlow"
	mv "${output_dir}/${2}/${1}/plasFlow/Unicycler_assemblies/${1}_uni_assembly/${1}_plasmid_assembly.fasta.TRIMMED.fasta" "${output_dir}/${2}/${1}/plasFlow/Unicycler_assemblies/${1}_uni_assembly/${1}_plasmid_assembly_trimmed.fasta"
	rm "${output_dir}/${2}/${1}/plasFlow/Unicycler_assemblies/${1}_uni_assembly/${1}_plasmid_assembly.fasta"
else
	echo "${output_dir}/${2}/${1}/Assembly/${1}_scaffolds_trimmed.fasta (Assembly) not found, cant do anything"
fi

#module unload PlasFlow/1.1
#module unload Python3/3.5.4
#module unload bowtie2/2.2.9
#module unload samtools/1.4.1
#module unload bam2fastq/1.1.0
#module unload Unicycler/0.4.4;
#module unload gcc/5.5;
#module unload SPAdes/3.11.1;
#module unload racon/1.2.0;

ml -Python3/3.5.4 -bowtie2/2.2.9 -samtools/1.4.1 -bam2fastq/1.1.0 -Unicycler/0.4.4 -SPAdes/3.13.0 -racon/1.3.1
