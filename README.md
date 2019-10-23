# Quaisar_singularity

Version of quaisar to be made into container

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
