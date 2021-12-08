# Quaisar_singularity
Quality, Assembly, Identification, Sequence type, Annotation, Resistance mechanisms for Hospital acquired infections (QuAISAR-H) is a mash-up of many publicly available tools with a splash of custom scripts with the purpose of producing a multi-layered quality checked report that identifies the taxonomy of and the Anti-microbial Resistance (AMR) elements from a paired end sequenced bacterial isolate.
This version uses containers to ease the necessity of having many preinstalled tools.

## Installation
    Dependencies - cURL, Singularity 3.7.2 and Go.
      Installation instructions for Singularity and go can be found at https://sylabs.io/guides/3.7/user-guide/quick_start.html#quick-installation-steps or https://github.com/jiangweiyao/Installing_Singularity#installing_singularity
      cURL can be installed by 'sudo apt-get install curl' or 'sudo yum install curl'

The script will install miniconda, if there is no version of conda already installed. It will then install an environment that contains Python3.6. When the pipeline is run it will activate the environment and will deactivate it when complete.

    To install
        A. Clone repository
        B. Navigate into the installation folder of the repository
        C. Run the instllation.sh script with the following parameters
            1.(-i) location of where to put the scripts to run the pipeline (installation directory)
            2.(-d) location of where to put the databases needed to run the pipeline (3.6 Gbs)
            3.(-w) location of where you expect to keep the output of the pipeline. This will be a folder that houses all runs (each within its own folder)

## Configuration
    There is a config.sh file within the script folder that has many configurable options
    The first section is system specific that gives the locations of where to find the scripts, databases, and output folders. These were set up during installation, but can be altered if desired. The installation program also attempted to learn the available number of processors to make sure the pipeline can run as quickly as possible. There are a large number of run-time options used by many tools that are defined within this file. The values of these options are the defaults used by the authors, however they may be adjusted to suit the problem.

## Running the pipeline
    This script will run each isolate through the pipeline in a serial fashion. The QuAISAR-H pipeline takes a folder of paired end reads and outputs multiple reports to the designated folder.

    Reads filenames need to have a postfix like one of these: _S*_L001_R*_00*.fastq[.gz], _S*_R*_0*X.fastq[.gz], _RX_00*.fastq[.gz], _[R]*.fastq[.gz]

    To run the pipeline use the following command with these parameters:
        A. ./quaisar_singularity.sh
            1. -i full path to the folder of paired-end reads
            2. -p run/set/project name for the set of reads being analyzed
            3. (optional) -o output_path_where_to_put_the_run/set/project_folder (if different than what was set during installation with the -w flag)
    Example: ./quaisar_singularity.sh -i /path/to/reads/folder -p run/set/project_name        

## Output
### Each run of the pipeline will produce the following files in the main (run/set/project name) folder
        1. A folder for each isolate's output files
        2. .log - standard out and err of all tools are directed to this file as well as being shown on the terminal
        3. _command.log - shows all singularity commands that were called during the run (and what the parameters were)
        4. .sum - is a concatenated version of all the individual isolates pipeline_stats files
        5. Seqlog_output.txt - A file that the authors use for run quality record keeping
        6. config_*.sh - copy of the config file used during the run.

### Each isolate folder will contain the following files and folders
    All Isolates -
    1. 16s folder
    2. ANI folder
    3. Assembly folder
    4. Assembly_Stats folder
    5. BUSCO folder
    6. c-sstar folder
    7. FASTQs folder
    8. GAMA folder
    9. gottcha folder
    10. kraken folder
    11. MLST folder
    12. plasmidFinder folder
    13. preQCcounts folder
    14. prokka folder
    15. removedAdapters folder
    16. srst2 folder
    17. trimmed
    18. _pipeline_stats.txt - summary and status of all steps performed in the pipeline
    19. .tax - determined taxonomy of isolate
    20. _time_summary.txt - estimate of length to complete each task



## Table with all external tools and versions used, along with example commands for each

|	Tool	|	Function	|	Version	|	command	|	command 2	|	Notes	|
|	---	|	---	|	---	|	---	|	---	|	---	|
|	BBDuk	|	Remove PhiX reads	|	BBMap(38.94)	|	bbduk.sh - Xmx20g threads=4 in=raw_R1.fastq in2=raw_R2.fastq out=noPhiX_R1.fsq out2=noPhiX_R2.fsq ref=phiX_adapter.fasta k=31 hdist=1	|		|		|
|	FastP	|	Remove illumina adapters and filter by quality	|	0.23.1	|	fastp -w 4 -i trimmed-noPhiX-R1.fsq -I trimmed-noPhiX-R2.fsq -o trimmed_R1_001.paired.fq --unpaired1 trimmed.single1.fq -O trimmed_R2_001.paired.fq --unpaired2 trimmed.single2.fq --adapter_fasta adapters.fasta -r --cut_right_window_size 20 --cut_right_mean_quality 30 -l 50 -g -5 20 -3 20 SLIDINGWINDOW:20:30 LEADING:20 TRAILING:20 MINLEN:50	|		|		|
|	Kraken	|	Taxonomic Identification/Contamination Detection	|	1.1.1	|	Reads: kraken --paired --db kraken_mini_db_location --preload --fastq-input --threads 4 --output sample_name.kraken --classified-out sample_name.classified  trimmed_R1_001.paired.fq trimmed_R2_001.paired.fq	|	Assembly: kraken --db kraken_mini_db_location --preload --threads 4 --output sample_name.kraken --classified-out sample_name.classified trimmed_assembly.fasta	|		|
|	SPAdes	|	Assembly	|	3.15.3	|	spades.py --careful --memory 32 --only-assembler --pe1-1 trimmed_R1_001.paired.fq --pe1-2 trimmed_R2_001.paired.fq --pe1-s trimmed.single.fq" -o output_directory --phred-offset 33 -t 12	|		|		|
|	QUAST	|	Assembly Quality	|	5.0.2	|	Quast.py -o output_directory trimmed_assembly.fasta	|		|		|
|	Prokka	|	Annotation	|	1.14.6	|	prokka --outdir output_directory trimmed_assembly.fasta	|		|		|
|	BUSCO	|	Determine quality of assembly and identification	|	3.0.2	|	run_BUSCO.py -i prokka_output_directory/sample_name.faa -o sample_name -l location_of_database -m prot	|		|	· Proper database is determined by matching lowest matching taxonomy to available databases	|
|	pyANI	|	Taxonomic Identification	|	0.2.11	|	average_nucleotide_identity.py -i directory_of_fastas -o output_directory --write_excel	|		|	· Directory of fastas contain the 20 closest genera matches based on mashtree distances	|
|	c-SSTAR	|	Anti-microbial Resistance Mechanism identification on Assembly	|	1.1.01	|	python3 c-SSTAR_gapped.py -g trimmed_assembly.fasta -s 98 -d  AR_database_location > sample_name.gapped_98.sstar	|		|		|
|	GAMMA	|	Anti-microbial Resistance Mechanism identification on Assembly	|	1.4	|	python3 GAMMA.py trimmed_assembly.fasta AR_database_location output_gama	|		|		|
|	SRST2	|	Anti-microbial Resistance Mechanism Identification on reads, Sequence Typing	|	0.2.0	|	AR: SRST2--input_pe trimmed_R1.fastq.gz trimmed_R2_001.fastq.gz --output output_directory –threads 4 --gene_db AR_datbase_location	|	MLST: SRST2--input_pe trimmed_R1.fastq.gz trimmed_R2_001.fastq.gz --output output_directory –threads 4 --mlst_db location_of_mlst_database --mlst_definitions location_of_MLST_definitions --mlst_delimiter MLST_definitions_file_delimiter	|	· Newest MLST database and definitions are downloaded as part of the script. The mlst delimiter is determined  using an included script within the SRST2, getmlst,  that must be run prior to SRST2	|
|	MLST	|	Sequence Typing	|	2.19.0	|	mlst trimmed_assembly.fasta > sample_name.mlst	|	mlst --scheme database_name trimmed_assembly.fasta > sample_name_database_name.mlst	|		|
|	Barrnap	|	Taxonomic Identification	|	0.9	|	barrnap --kingdom bac --threads 4  trimmed_assembly.fasta > rRNA_seqs.fasta	|		|		|
|	plasmidFinder	|	Anti-microbial Resistance Mechanism Identification on plasmid replicons	|	2.1.1	|	plasmidfinder -i trimmed_assembly.fasta -o output_directory -k 95.00 -p enterobacteriaceae|gram_positive	|		|		|
|   bowtie2 |   read aligner    |   2.2.9   |   bowtie2-build -f plasFlow_results.tsv_chromosomes.fasta bowtie2_sample_name_chr |   bowtie2 -x sample_name_chr -1 R1_001.paired.fq -2 R2_001.paired.fq -S sample_name.sam -p 12 --local |       |       |
|   samtools    |   sam converter   |   1.14    |   samtools view -bS sample_name.sam > sample_name.bam |   sort -n sample_name.bam -o sample_name.bam.sorted
|   bedtools    |   bam converter   |   2.30.0  |   bamToFastq -i sample_name.bam.sorted -fq sample_name_R1_bacterial.fastq -fq2 sample_name__R2_bacterial.fastq



##Flag table of output summaries

Title|ALERT|WARNING|FAILED
---|---|---|---
Time|-time_summary.txt missing||
FASTQs[R1-R2]||-Only one read file found|-No read files found or empty
QC counts||-<1,000,000 reads|-No sample_counts.txt found
Q30_R1%||-<90%|-No sample_counts.txt found
Q30_R2%||-<70%|-No sample_counts.txt found
BBDUK=PhiX[R1-R2]|-No reads removed|-More phix reads found than raw reads or no removedAdapters folder found
Trimming||-Only one trimmed read file found|No trimmed reads files found
QC count after trim||-<500,000 reads|-No trimmed_sample_counts.txt found
kraken preassembly|||-.kraken(.gz) missing
krona-kraken-preasmb|||-.krona or .html missing
Pre classfify||-unclassified reads >30%|-no classified reads or kraken_summary_paired.txt missing
pre Class Contam.|-More than one species found above 25% threshold|-No species found above 25%
Assembly|||-scaffolds.fasta is missing
Contig Trim||->200 contigs remain|-scaffolds_trimmed.fasta missing
kraken postassembly|||-.kraken(.gz) missing
krona-kraken-pstasmb|||-.krona or .html missing
post Classify||-unclassified reads >30%  or no single species >50%|-no classified reads or kraken_summary_assembled.txt missing
post Clas Contam.|-More than one species found above 25% threshold|-No species found above 25%
kraken weighted|||-.kraken missing
krona-kraken-weight|||-.krona or .html missing                                                         -kraken preassembly failed
weighted Classify|||-unclassified reads >30%|-no single species >50% or no classified reads or -kraken_summary_assembled_BP.txt missing
weighted Contam.|||-More than one species found above 25% threshold or no species found above 25%
QUAST|||report.tsv missing
Taxa|||No species was determined
Assembly ratio||-if taxonomy is not in DB|-ratio >1.2x or <.8x
Raw coverage|-Coverage >90x||-<40x
Trimmed coverage|-Coverage >90x||-<40x
prokka|||-prokka.gbf missing
BUSCO|||-<90% gene count match or -summary file not found
ANI_REFSEQ|REFSEQ database date is not current||-<95% match, <70% coverage,best_hits file or ANI folder missing
c-SSTAR|-NO known AMR genes present,database is not current||-summary.txt or c-sstar folder missing
GAMA|-NO known AMR genes present,database is not current||.GAMA or GAMA folder missing
SRST2|GAMA|-NO known AMR genes present,database is not current||genes file or srst2 folder missing
MLST|-No scheme found for taxa|-Type does not exist|-.mlst,1+ allele is missing,taxa not defined
MLST-srst2|-No scheme found for taxa, more than 2 srst2 files found, more than 1 srst2 file found for non A.Baum or E.coli|-Type does not exist|-.mlst,1+ allele is missing,taxa not defined
16s_best_hit||-species not found|-Genus not found, 16s_blast_id.txt missing,No reads found,Unclassifiable reads found
16s_largest_hit||-species not found|-Genus not found, 16s_blast_id.txt missing,No reads found,Unclassifiable reads found
plasmidFinder|||-results_table_summary.txt missing,plasmidFinder folder missing
