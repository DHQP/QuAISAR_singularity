#!/usr/bin/env python3

#
# Description: Conversion script from MetaPhlAn output to Krona text input file
#
# Usage: python ./project_parser.py -p metaphlan_input_file -k krona_output_file
#
# Output location: Parameter
#
# Modules required: None
#
# v1.0.2 (08/04/2021)
#
# Created by Nick Vlachos (nvx4@cdc.gov)
#

import sys
import glob
import math
import argparse

# Parse all arguments from command line
def parseArgs(args=None):
	parser = argparse.ArgumentParser(description='Script to consolidate all info into a single report for outbreak analysis')
	#parser.add_argument('-c', '--csstar', required=True, help='input summary filename')
	parser.add_argument('-p', '--plasmid', required=True, help='input plasmid filename')
	parser.add_argument('-s', '--summary', required=True, help='input sample summary file filename')
	parser.add_argument('-o', '--output', required=True, help='output csv filename')
	parser.add_argument('-d', '--database', required=True, help='database file used in AR discovery')
	parser.add_argument('-m', '--sim', required=True, help='database %id used for filtering')
	parser.add_argument('-l', '--length', required=True, help='database % length used for filtering')
	return parser.parse_args()


# main function that sorts and formats all AR genes found using csstar, GAMA and srst2 that have already been filtered for % identity and % length
def do_AR(input_summary_AR, input_plas, output_file, DB_name, percent_sim, percent_length):
	all_ARs_in_file=[]
	samples=[]
	input_summary=open(input_summary_AR,'r')
	summary_line = input_summary.readline().strip()
	counter=0
	#print("Start")
	while summary_line != '':
		#print(counter, summary_line)
		#print("Start csstar loop")
		summary_line_sections=summary_line.split("	")
		csstar_list=summary_line_sections[10].split(",")
		#print(csstar_list)
		ar_dict={}
		for csstar_gene in csstar_list:
			gene_name=csstar_gene.split("[")[0]
			if gene_name == "No AR genes discovered":
				gene_stats="[0/0:#-]C"
			elif gene_name == "NO CURRENT FILE":
				gene_name="NO-"+DB_name+"-CSSTAR-FILE at "+percent_sim+"% identity"
				gene_stats="[NA/NA:#-]C"
			else:
				gene_stats="["+csstar_gene.split("[")[1]+"C"
			ar_dict[gene_name]=gene_stats
			if gene_name not in all_ARs_in_file:
				all_ARs_in_file.append(gene_name)
				#print("Adding", gene_name)

		srst2_list=summary_line_sections[11].split(",")
		for srst2_gene in srst2_list:
			if srst2_gene == "No AR genes discovered":
				gene_name="No AR genes discovered"
				gene_stats="[0/0]S"
				#print("Looking up", gene_name, "in csstar dic")
				if ar_dict.get(gene_name):
					ar_dict[gene_name]=""+ar_dict.get(gene_name)+":"+gene_stats
				else:
					print("New AR-less isolate found in srst2")
					ar_dict[gene_name]=gene_stats
				if gene_name not in all_ARs_in_file:
					all_ARs_in_file.append(gene_name)
			elif srst2_gene == "NO CURRENT FILE":
				gene_name="NO-"+DB_name+"-SRST2-FILE"
				gene_stats="[NA/NA]S"
				#print("Looking up", gene_name, "in csstar dic")
				if ar_dict.get(gene_name):
					ar_dict[gene_name]=""+ar_dict.get(gene_name)+":"+gene_stats
				else:
					print("New AR-less isolate found in srst2")
					ar_dict[gene_name]=gene_stats
				if gene_name not in all_ARs_in_file:
					all_ARs_in_file.append(gene_name)
			else:
				gene_name=srst2_gene.split("[")[0]
				gene_stats="["+srst2_gene.split("[")[1]+"S"
				#print("Looking up", gene_name, "in csstar dic")
				if ar_dict.get(gene_name):
					if ar_dict.get(gene_name) != "No Other AR genes":
						#print("Found", gene_name, "in both outputs")
						#print("New value: "+ar_dict.get(gene_name)+":"+gene_stats)
						ar_dict[gene_name]=""+ar_dict.get(gene_name)+":"+gene_stats
					#else:
						#print("No AR found in csstar for", gene_name)
					#	:
				else:
					#print("New gene", gene_name,"found in srst2")
					ar_dict[gene_name]=gene_stats
				if gene_name not in all_ARs_in_file:
					all_ARs_in_file.append(gene_name)
					#break

		GAMA_list=summary_line_sections[12].split(",")
		for GAMA_gene in GAMA_list:
			if GAMA_gene == "No AR genes discovered":
				gene_name="No AR genes discovered"
				gene_stats="[0/0]G"
				#print("Looking up", gene_name, "in csstar dic")
				if ar_dict.get(gene_name):
					ar_dict[gene_name]=""+ar_dict.get(gene_name)+":"+gene_stats
				else:
					print("New AR-less isolate found in GAMA")
					ar_dict[gene_name]=gene_stats
				if gene_name not in all_ARs_in_file:
					all_ARs_in_file.append(gene_name)
			elif GAMA_gene == "NO CURRENT FILE":
				gene_name="NO-"+DB_name+"-GAMA-FILE"
				gene_stats="[NA/NA]G"
				#print("Looking up", gene_name, "in csstar dic")
				if ar_dict.get(gene_name):
					ar_dict[gene_name]=""+ar_dict.get(gene_name)+":"+gene_stats
				else:
					print("New AR-less isolate found in GAMA")
					ar_dict[gene_name]=gene_stats
				if gene_name not in all_ARs_in_file:
					all_ARs_in_file.append(gene_name)
			else:
				gene_name=GAMA_gene.split("[")[0]
				gene_stats="["+GAMA_gene.split("[")[1]+"G"
				#print("Looking up", gene_name, "in csstar dic")
				if ar_dict.get(gene_name):
					if ar_dict.get(gene_name) != "No Other AR genes":
						#print("Found", gene_name, "in both outputs")
						#print("New value: "+ar_dict.get(gene_name)+":"+gene_stats)
						ar_dict[gene_name]=""+ar_dict.get(gene_name)+":"+gene_stats
					#else:
						#print("No AR found in csstar for", gene_name)
					#	:
				else:
					#print("New gene", gene_name,"found in srst2")
					ar_dict[gene_name]=gene_stats
				if gene_name not in all_ARs_in_file:
					all_ARs_in_file.append(gene_name)
			#break
		samples.append([summary_line_sections[0], summary_line_sections[1], summary_line_sections[2], summary_line_sections[3],  summary_line_sections[4], summary_line_sections[5], summary_line_sections[6], summary_line_sections[7], summary_line_sections[8], summary_line_sections[9], ar_dict])
		#print("Total AR genes in sample set:", len(all_ARs_in_file)-1)
		summary_line = input_summary.readline().strip()
	input_summary.close
	all_ARs_in_file.sort()
	if len(all_ARs_in_file) == 0:
		print("\n")
		print("Total AR genes in sample set: 0")
	else:
		print("Total AR genes in sample set:",len(all_ARs_in_file))
		#print(*all_ARs_in_file, sep = "\n")
		for gene in all_ARs_in_file:
			if gene != "No other AR genes":
				print (gene)
	print()

	# print("1#", len(samples))
	# for sample in samples:
	# 	print(sample[1])

	#Parse plasmid summary file
	all_plasmids_in_file=[]
	plas_file=open(input_plas, 'r')
	line = plas_file.readline().strip()
	sample_p_plasmids_dict={}
	sample_f_plasmids_dict={}
	current_id=""
	counter=0
	while line != '':
		#print(counter, line)
		plasmid_line_sections=line.split("	")
		#print("Current id:", current_id, ":", plasmid_line_sections[0], plasmid_line_sections[1])

		if plasmid_line_sections[0]+"/"+plasmid_line_sections[1] != current_id and current_id != "":
			#print("New name! - ")
			for sample_index in range(0,len(samples)):
				#print("Looking for", current_id, ", found", samples[sample_index][0].strip()+"/"+samples[sample_index][1].strip())
				if current_id == samples[sample_index][0]+"/"+samples[sample_index][1].strip():
					# if current_id == "MUGSI_197/MUGSI-99_FAILED":
					# 	print(samples[sample_index], "adding full plasmids")
					# 	fcount = 0
					# 	pcount = 0
					# 	for key, value in sample_f_plasmids_dict.items():
					# 		print (key, value)
					# 		fcount += 1
					# 	if fcount == 0:
					# 		print("adding nothing {}")
					# 	samples[sample_index].append(sample_f_plasmids_dict)
					# 	print(samples[sample_index], "adding plasmid plasmids")
					# 	for key, value in sample_p_plasmids_dict.items():
					# 		print (key, value)
					# 		pcount += 1
					# 	if pcount == 0:
					# 		print("adding nothing {}")
					# 	samples[sample_index].append(sample_p_plasmids_dict)
					# 	print(len(samples[sample_index]), samples[sample_index])
					# else:
					samples[sample_index].append(sample_f_plasmids_dict)
					samples[sample_index].append(sample_p_plasmids_dict)
					break
				#else:
					#print("Sample", current_id, "does not exist")
			current_id=plasmid_line_sections[0]+"/"+plasmid_line_sections[1].strip()
			#print(current_id)
			sample_f_plasmids_dict={}
			sample_p_plasmids_dict={}


		source_assembly=plasmid_line_sections[2]
		#print("Test:"+plasmid_line_sections[4]+":")
		#if plasmid_line_sections[4].find("_contigs-") >= 0:
		#	line = plas_file.readline().strip()
		#	continue

		if current_id == "":
			current_id = plasmid_line_sections[0]+"/"+plasmid_line_sections[1]

		if plasmid_line_sections[3] == "No_Plasmids_Found":
			plas_perc_id="-"
			plas_perc_length="-"
		else:
			plas_perc_id=math.floor(float(plasmid_line_sections[4]))
			#print("testing:", plasmid_line_sections[5].split("/")[0], plasmid_line_sections[5].split("/")[1])
			plas_perc_length=(100*int(plasmid_line_sections[5].split("/")[1])//int(plasmid_line_sections[5].split("/")[0]))
		plas_match_info="["+str(plas_perc_id)+"/"+str(plas_perc_length)+"]"
		if source_assembly == "full_assembly":
			#print("Adding:", plasmid_line_sections[3], "to sample_f_plasmids of", current_id)
			sample_f_plasmids_dict[plasmid_line_sections[3]]="["+str(plas_perc_id)+"/"+str(plas_perc_length)+"]"
		elif source_assembly == "plasmid_assembly":
			#print("Adding:", plasmid_line_sections[3], "to sample_p_plasmids of", current_id)
			sample_p_plasmids_dict[plasmid_line_sections[3]]="["+str(plas_perc_id)+"/"+str(plas_perc_length)+"]"
		if len(plasmid_line_sections) > 1:
			if plasmid_line_sections[3] not in all_plasmids_in_file:
				all_plasmids_in_file.append(plasmid_line_sections[3])
				#print("Adding to project list", plasmid_line_sections[3])
		#else:
			#print("Line length:", len(csstar_plasmid_line_sections))
		#print()
		line = plas_file.readline().strip()
		counter=counter+1

	for sample_index in range(0,len(samples)):
		#print("Looking for", current_id, ", found", samples[sample_index][0].strip()+"/"+samples[sample_index][1].strip())
		if current_id == samples[sample_index][0]+"/"+samples[sample_index][1].strip():
			#print(samples[sample_index], "adding", sample_f_plasmids_dict, "and", sample_p_plasmids_dict)
			samples[sample_index].append(sample_f_plasmids_dict)
			samples[sample_index].append(sample_p_plasmids_dict)
			break
		else:
			print("Sample", current_id, "does not exist")
	all_plasmids_in_file.sort()
	if len(all_plasmids_in_file) == 0:
		#print("\n")
		print("Total plasmid replicons in sample set: 0")
		#all_plasmids_in_file=["No_Plasmids_Found"]
	elif len(all_plasmids_in_file) == 1 and all_plasmids_in_file[0] == "No_Plasmids_Found":
		print("Total plasmid replicons in sample set: 0")
	else:
		print("Total plasmid replicons in sample set:", len(all_plasmids_in_file))
		print(*all_plasmids_in_file, sep= "\n")
	print()
	plas_file.close

	all_ar_and_plasmids=all_ARs_in_file+["|"]+all_plasmids_in_file
	#all_AR_to_write=all_ARs_in_file
	#all_AR_to_write.insert(0,",")
	#all_AR_to_write.insert(0,",")
	#all_AR_to_write=','.join(map(str, all_AR_to_write))
	header="id,Path__autocolour,Project__autocolour,Species__autocolour,Species_determinant__autocolour,Species_Support__autocolour,MLST_Pasteur__autocolour,MLST_Pasteur_alleles__autocolour,ALT_MLST__autocolour,ALT_MLST_alleles__autocolour,AR_Database__autocolour"
	for thing in all_ar_and_plasmids:
		header = header + "," + thing + "__autocolour"
	#header = header[:-1]
	summary_out=open(output_file, 'w')
	summary_out.write(header+'\n')
	#all_AR_to_write=all_AR_to_write[2:]
	#print("List:", all_ar_and_plasmids)
	#for sample in samples:
	#	print ("2:",sample[0])
	#return
	for sample in samples:
		sample_details=[sample[1], sample[2], sample[0], sample[3], sample[4], sample[5], sample[6], sample[7], sample[8], sample[9], DB_name+"(["+percent_sim+"NT/"+percent_length+"]C:[XNT/"+percent_sim+"AA/"+percent_length+"]G:["+percent_sim+"NT/"+percent_length+"]S)"]
		#print("pre:", len(sample),sample)
		for gene in all_ar_and_plasmids:
			status=" "
			#print(sample_details[2], sample_details[3], gene, len(sample))
			#if len(sample_details) == 11:
			#	print(sample)
			if gene == "|":
				sample_details.append(gene)
				continue
			if sample[10].get(gene):
				status=sample[10].get(gene)
			elif sample[11].get(gene):
				if sample[12].get(gene):
					status="F:"+sample[11].get(gene)+";P:"+sample[12].get(gene)
				else:
					status="F:"+sample[11].get(gene)
			elif sample[12].get(gene):
				status="P:"+sample[12].get(gene)
			sample_details.append(status)
		#print("Post Sample check", sample_details)
		sample_details=','.join(map(str, sample_details))
		summary_out.write(sample_details+"\n")
	summary_out.close
#print (sys.argv[1:])



print("Parsing project AR files ...\n")
args = parseArgs()
do_AR(args.summary, args.plasmid, args.output, args.database, args.sim, args.length)
