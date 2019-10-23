# Quaisar_singularity

Version of quaisar to be made into container

Table with all Tools and versions used, along with example commands for each
Tool	|	Function/command	|	Version	|	Newest  	|	command	|	command 2	|	Notes
---	|	---	|	---	|	---	|	---	|	---	|	---
BBDuk 	|	Remove PhiX reads	|	BBMap(38.26)	|	38.26(38.42)	|	bbduk.sh - Xmx20g threads=12 in=raw_R1.fastq in2=raw_R2.fastq out=noPhiX_R1.fsq out2=noPhiX_R2.fsq ref=phiX_adapter.fasta k=31 hdist=1	|		|	
Trimmomatic	|	Remove illumina adapters and filter by quality	|	0.35 	|	0.36 (0.38)	|	trimmomatic PE -phred33 -threads 12 noPhiX_R1.fsq noPhiX_R2.fsq trimmed_R1_001.paired.fq trimmed_R1_001.unpaired.fq trimmed_R2_001.paired.fq trimmed_R2_001.unpaired.fq ILLUMINACLIP:adapters.fasat:2:30:10:8:TRUE SLIDINGWINDOW:20:30 LEADING:20 TRAILING:20 MINLEN:50	|		|	
Kraken	|	Taxonomic Identification/Contamination Detection	|	0.10.5	|	2.0.8(2.0.7)	|	Reads: kraken --paired --db kraken_mini_db_location --preload --fastq-input --threads 12 --output sample_name.kraken --classified-out sample_name.classified  trimmed_R1_001.paired.fq trimmed_R2_001.paired.fq	|	Assembly: kraken --db kraken_mini_db_location --preload --threads 14 --output sample_name.kraken --classified-out sample_name.classified trimmed_assembly.fasta	|	
Gottcha	|	Taxonomic Identification (Species database)	|	1.0b	|	1.0b	|	gottcha.pl --mode all --outdir output_directory --input paired.fq --database location_of_gottcha_database	|		|	· Paired.fq is the concatenated file of trimmed R1 and R2 read files
SPAdes	|	Assembly	|	3.13.0	|	3.13.0	|	spades.py --careful --memory 32 --only-assembler --pe1-1 trimmed_R1_001.paired.fq --pe1-2 trimmed_R2_001.paired.fq --pe1-s trimmed.single.fq" -o output_directory --phred-offset 33 -t 12	|		|	
QUAST	|	Assembly Quality	|	4.3	|	5.0(5.0.2)	|	Quast.py -o output_directory trimmed_assembly.fasta	|		|	
Prokka	|	Annotation	|	1.12	|	1.13.3	|	prokka --outdir output_directory trimmed_assembly.fasta	|		|	
BUSCO	|	Determine quality of assembly and identification	|	3.0.1	|	3.0.1(3.1.0)	|	run_BUSCO.py -i prokka_output_directory/sample_name.faa -o sample_name -l location_of_database -m prot	|		|	· Proper database is determined by matching lowest matching taxonomy to available databases
pyANI	|	Taxonomic Identification	|	0.2.7	|	0.2.7(0.2.9)	|	average_nucleotide_identity.py -i directory_of_fastas -o output_directory --write_excel	|		|	· Directory of fastas contain the 20 closest genera matches based on mashtree distances
c-SSTAR	|	Anti-microbial Resistance Mechanism identification on Assembly	|	1.1.01	|	1.1.01	|	Normal: python3 c-SSTAR_gapped.py -g trimmed_assembly.fasta -s 98-d  AR_database_location > sample_name.gapped_98.sstar	|	Plasmid: python3 c-SSTAR_gapped.py -g plasmid_assembly.fasta -s 40-d  AR_database_location > sample_name.gapped_40.sstar	|	
SRST2	|	Anti-microbial Resistance Mechanism Identification on reads, Sequence Typing	|	0.2.0	|	0.2.0	|	AR: SRST2--input_pe trimmed_R1.fastq.gz trimmed_R2_001.fastq.gz --output output_directory –threads 12 --gene_db AR_datbase_location	|	MLST: SRST2--input_pe trimmed_R1.fastq.gz trimmed_R2_001.fastq.gz --output output_directory –threads 12 --mlst_db location_of_mlst_database --mlst_definitions location_of_MLST_definitions --mlst_delimiter MLST_definitions_file_delimiter	|	· Newest MLST database and definitions are downloaded as part of the script. The mlst delimiter is determined  using an included script within the SRST2, getmlst,  that must be run prior to SRST2
MLST	|	Sequence Typing	|	2.16	|	2.16	|	mlst trimmed_assembly.fasta > sample_name.mlst	|	mlst --scheme  database_name trimmed_assembly.fasta > sample_name_database_name.mlst	|	
Barrnap	|	Taxonomic Identification	|	0.8	|	0.8	|	barrnap --kingdom bac --threads 12  trimmed_assembly.fasta > rRNA_seqs.fasta	|		|	
plasmidFinder	|	Anti-microbial Resistance Mechanism Identification on plasmid replicons	|	1.3	|	1.3(2.0)	|	plasmidfinder -i trimmed_assembly.fasta -o output_directory -k 95.00 -p enterobacteriaceae|gram_positive	|		|	
FASTQ_quality_printer.py	|	Quality Control for Reads	|	N/A	|	N/A	|	Initial read quality: python2 Fastq_Quality_Printer.py -1 R1.fastq -2 R2.fastq > counts.txt	|	trimmed read quality: python2 Fastq_Quality_Printer.py -1 trimmed_R1.fastq -2 trimmed_R2.fastq > trimmed_counts.txt	|	
removeShortContigs.py	|	Quality Control for Assembly	|	N/A	|	N/A	|	Assembly trimming: python3 removeShortContigs.py -i trimmed_assembly.fasta -t 500 -s normal_SPAdes	|	Plasmid assembly trimming: python3 removeShortContigs.py -i trimmed_plasmid_assembly.fasta -t 2000 -s normal_SPAdes	|	


Flag table of output summaries
Title|ALERT|WARNING|FAILED
---|---|---|---
Time|-time_summary.txt missing ||
FASTQs||-Only one read file found|-No reads files found
QC counts||-<1,000,000 reads|-No sample_counts.txt found
Q30_R1%||-<90%|-No sample_counts.txt found
Q30_R2%||-<70%|-No sample_counts.txt found
Trimming||Only one trimmed read file found|No trimmed reads files found
QC count after trim||-<500,000 reads|-No trimmed_sample_counts.txt found
kraken preassembly|||-.kraken missing
krona-kraken-preasmb|||-.krona or .html missing                                                        -kraken preassembly failed
Pre classfify||-unclassified reads >30%|-no classified reads                                               -kraken_summary_paired.txt missing
GOTTCHA_S||-.tsv OR .html missing|-Both .tsv and .html missing
GottchaV1 Classifier||-unclassified reads >30%|-no classified reads                                              -gottcha_species_summary.txt missing
Assembly|||-scaffolds.fasta is missing
plasmid Assembly|||-Plasmid folder is missing
Contig Trim||->200 contigs remain|-scaffolds_trimmed.fasta missing
Plasmids contig Trim|||-plasmid_scaffolds_trimmed.fasta
kraken postassembly|||-.kraken missing
krona-kraken-pstasmb|||-.krona or .html missing                                                                      -kraken preassembly failed
post Classify||-unclassified reads >30%                     -no single species >50%|-no classified reads                                                                   -kraken_summary_assembled.txt missing
kraken weighted|||-.kraken missing
krona-kraken-weight|||-.krona or .html missing                                                         -kraken preassembly failed
weighted Classify||-unclassified reads >30%                    -no single species >50%|-no classified reads                                                                                        -kraken_summary_assembled_BP_data.txt missing
QUAST|||report.tsv missing
Assembly ratio||-if taxonomy is not in DB|-ratio >1.2x or <.8x
Raw coverage|-Coverage >90x||-<40x
Trimmed coverage|-Coverage >90x||-<40x
prokka|||-prokka.gbf missing
BUSCO|||-<90%                                                                                            -summary file not found
ANI|||-<95%                                                                                         -best_hits or ANI/ folder missing
c-SSTAR|-NO known AMR genes present                        -ran against older DB||-summary.txt missing
c-SSTAR_plasmid|-NO known AMR genes present                        -ran against older DB||-summary.txt missing
SRST2|-NO known AMR genes present                       -ran against older DB||summary.txt missing
MLST||-Type does not exist|-file or scheme does not exist
16s_best_hit||-species not found|-Genus not found, 16s_blast_id.txt missing
16s_largest_hit||-species not found|-Genus not found, 16s_blast_id.txt missing
plasmid|||-results_table_summary.txt missing                                              -plasmid/ missing
plasmid-plasmidAsmb|||-results_table_summary.txt missing                                    -plasmid_on_plasmidAssembly/ missing

