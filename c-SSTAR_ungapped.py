#!/usr/bin/env python

#$ -o C-SSTAR_u.out
#$ -e C-SSTAR_u.err
#$ -N C-SSTAR_u
#$ -cwd
#$ -q short.q

#
# Description: Script to discover best AR gene hits from each cluster from an srst2 formatted database not allowing gaps in matches
#
# Usage is python ./c-SSTAR_gapped.py -d srst2_formatted_database -g fasta_genome -b basename_of_output_file -o output_directory -s similarity_%_minimum [-e number_of_bases_at _contig_edge] [-v show_version]
#
# Output location: parameter
#
# Modules required: Biopython must be available in python instance
#
# Created by Chris Gulvik, Tom de Man and Adrian Lawsin (kqj9@cdc.gov)
#

__version__ = '1.2c'

import argparse
import logging
import os
import pwd
import re
import subprocess
import sys
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

def parseArgs(args=None):
	parser = argparse.ArgumentParser(description='c-SSTAR is a CLI utility for rapidly identifying antibiotic resistance gene determinants in bacterial genomes')
	parser.add_argument('-d', '--database', required=True, help='a SSTAR-formatted FastA database of AR gene sequences')
	parser.add_argument('-g', '--genome', required=True, help='a FastA genome')
	parser.add_argument('-b', '--basename', help='basename of output files')
	parser.add_argument('-o', '--outdir', default=os.getcwd(), help='output directory [default: cwd]')
	parser.add_argument('-s', '--similarity', type=int, default=95, help='minimum percent nucleotide similarity [default: 95]')
	parser.add_argument('-e', '--edge', type=int, default=50, help='number of bases at each contig edge to report as \'$\' for end [default 50]')
	parser.add_argument('-v', '--version', action='version', version='%(prog)s v{}'.format(__version__))
	return parser.parse_args()

def syscall(syscmd):
	with open(os.devnull) as dump:
		returncode = subprocess.call(syscmd, stdout=dump, stderr=dump, shell=True)
		if returncode != 0:
			logging.error('failed syscall ' + syscmd)
			sys.exit('ERROR: failed syscall ' + syscmd)

def translateSeq(nuclSeq):
	'''use biopython to take in a nucleotide sequence (string) and return a protein sequence'''
	proteinSeq = Seq(nuclSeq, generic_dna).translate(cds=False, to_stop=False, stop_symbol='*')
	return proteinSeq

def internalSTOPcodon(candidate):
	'''counts number of internal stop codons; requires sequence from database to be in frame'''
	if '-' in candidate[11]:
		nucSeq = candidate[11].replace('-', '')
	else:
		nucSeq = candidate[11]
	if int(int(candidate[8])%3) == 1:
		frameStart = int(0)
	elif int(int(candidate[8])%3) == 0:
		frameStart = int(1)
	elif int(int(candidate[8])%3) == 2:
		frameStart = int(2)
	if frameStart > 0:
		nucSeq = nucSeq[frameStart:]
	if len(nucSeq)%3 == 0:
		protein = translateSeq(nucSeq)
	elif len(nucSeq)%3 > 0:
		protein = translateSeq(nucSeq[:-(len(nucSeq)%3)])
	else:
		sys.exit('ERROR: incorrect nucleotide length ({}) after trimming'.format(len(nucSeq)))
	numInternalSTOPcodons = len(re.findall(r'\*[ABCDEFGHIKLMNPQRSTVWYZ]', str(protein)))
	return (protein, str(numInternalSTOPcodons))

def tagHit(l, edge):
	if int(l[5]) != int(l[7]):  #incomplete len; SRST2='? indicates that there was uncertainty in at least one of the alleles'
		l = [ l[0], l[1]+'?', l[2]+'?', l[3], l[4], l[5], l[6], l[7], l[8], l[9], l[10], l[11], l[12], l[13], l[14], l[15] ]
	if int(l[5]) == int(l[7]) and int(l[4]) < 100:  #full length match but pident!=100%; SRST2='* [...] indicates that there were mismatches against at least one of the alleles. This suggests that you have a novel variant [...] rather than a precise match'
		l = [ l[0], l[1]+'*', l[2]+'*', l[3], l[4], l[5], l[6], l[7], l[8], l[9], l[10], l[11], l[12], l[13], l[14], l[15]]
	if (int(l[8]) - edge) < 0:  #Test for edge hits on left edge
		if (int(l[7]) - int(l[9])) > 0:  #require incomplete alignlen at edge for '$' designation (distinguishing ^ from $ is not necessary so KISS)
			l = [ l[0], l[1]+'$', l[2]+'$', l[3], l[4], l[5], l[6], l[7], l[8], l[9], l[10], l[11], l[12], l[13], l[14], l[15] ]
	if (int(l[9]) + edge) > int(l[10]):  #Test for edge hits on right edge
		if (int(l[7]) + int(l[8])) > int(l[10]):  #require incomplete alignlen
			l = [ l[0], l[1]+'$', l[2]+'$', l[3], l[4], l[5], l[6], l[7], l[8], l[9], l[10], l[11], l[12], l[13], l[14], l[15] ]
	(prot, numSTOP) = internalSTOPcodon(l)
	if int(numSTOP) > 0:  #TR indicates 'truncated protein translation'
		l = [ l[0], l[1]+'-TRUNC', l[2]+'-TRUNC', l[3], l[4], l[5], l[6], l[7], l[8], l[9], l[10], l[11], l[12], l[13], l[14], l[15] ]
	return l

def main(args=None):
	args = parseArgs()
	outdir = args.outdir
	genome = args.genome
	if args.basename is not None:
		baseGenome = args.basename
	else:
		baseGenome = os.path.splitext(os.path.basename(genome))[0]
	if not os.path.exists(outdir):
		os.mkdir(outdir)
	logging.basicConfig(filename='{}/c-SSTAR_{}.log'.format(outdir, baseGenome), format='%(asctime)s: %(levelname)s: %(message)s', datefmt='%d-%m-%Y %I:%M:%S %p', level=logging.INFO)
	logging.info('c-SSTAR version: {}'.format(__version__))
	logging.info('user: {}'.format(pwd.getpwuid(os.getuid()).pw_name))
	logging.info('release: {}'.format(os.uname()[3]))
	logging.info('shell env: {}'.format(pwd.getpwuid(os.getuid()).pw_shell))
	logging.info('cwd: {}'.format(pwd.getpwuid(os.getuid()).pw_dir))
	logging.info('python version: {}'.format(sys.version))
	logging.info(subprocess.check_output('command -v blastn', shell=True).rstrip())
	logging.info(subprocess.check_output('blastn -version | tail -n1', shell=True).rstrip())
	database = args.database
	similarity = args.similarity
#	print("in:",genome,"out-",os.path.join(outdir, baseGenome))
	syscall('makeblastdb -in {} -out {} -dbtype nucl'.format(genome, os.path.join(outdir, baseGenome)))
	syscall('blastn -task blastn -query {0} -db {out} -out {out}.blastn.tsv -ungapped -evalue 1e-5 -max_target_seqs 1 -perc_identity {1} -culling_limit 1 -outfmt "6 qseqid sseqid pident length mismatch gaps qstart qend sstart send evalue bitscore qlen slen sseq"'.format(database, similarity, out=os.path.join(outdir, baseGenome)))
	os.remove(os.path.join(outdir, baseGenome + '.nin'))
	os.remove(os.path.join(outdir, baseGenome + '.nsq'))
	os.remove(os.path.join(outdir, baseGenome + '.nhr'))
	#print("test")
	#print (""+str(sys.version_info[0])+"")
	#print(os.path.join(outdir, baseGenome + '.blastn.tsv'))
	with open(os.path.join(outdir, baseGenome + '.blastn.tsv')) as infile:
		best = ['-1','a','a','a',0,'-1','a',0,'a','a',0, 'a']
		currentClusterNr = '-1'
		topHits = []
		for l in infile:
			blastOut = [x for x in l.split('\t')]
			# Special for cdiff project
			thresHold = (int(blastOut[12])/5)
			# Original
#			thresHold = (int(blastOut[12])/5)*2
			if int(blastOut[3]) > thresHold:  #>40% overlap required
				enzymeParts = [s for s in blastOut[0].split('__')]
				clusterNr = enzymeParts[0]
				ident, dec = blastOut[2].split('.')
				pident = int(ident)
				bitscore = blastOut[11]
				#print("blastOut is", len(blastOut), "elements")
				#print (", ".join(map(str, blastOut)))
				#print("enzymeParts is", len(enzymeParts), "elements")
				#print (", ".join(map(str, enzymeParts)))
				#print()
				#for i in range(0, len(blastOut)):
				#	print("B:", i, blastOut[i])
				#for i in range(0, len(enzymeParts)):
				#	print("E:", i, enzymeParts[i])
				#print("\n\n\n\n\n")
				#candidate = [clusterNr, enzymeParts[1], enzymeParts[2], blastOut[1] , pident, blastOut[3], bitscore, blastOut[12], blastOut[6], blastOut[7], blastOut[13], blastOut[14].rstrip(), enzymeParts[4] , enzymeParts[5], round(100*(int(blastOut[12])/int(blastOut[3]))), blastOut[4], ]
				# For use with reporting SNPs
				candidate = [clusterNr, enzymeParts[1], enzymeParts[2], blastOut[1] , pident, blastOut[3], bitscore, blastOut[12], blastOut[6], blastOut[7], blastOut[13], blastOut[14].rstrip(), enzymeParts[4] , enzymeParts[5], int(100*int(blastOut[3])/int(blastOut[12])), blastOut[4]]

				if clusterNr == currentClusterNr:
					if float(best[6]) < float(candidate[6]):
						best = candidate
				else:
					if best[4] >= similarity:
						topHits.append(best)
					currentClusterNr = clusterNr
					best = candidate
		if best[4] >= similarity:
			topHits.append(best)
	for i in topHits:
		hit = tagHit(i, args.edge)
		# Reports numbers of mismatches at the end
		#print (str(hit[13]) + '\t' + str(hit[12]) + '\t' + hit[1] + '\t' + hit[2] + '\t' + hit[3] + '\t' + str(hit[4]) + '%\t' + str(hit[5]) + '\t' + str(hit[7]) + '\t' + str(hit[14]) + '\t' + str(hit[15]))
		# Does not report SNPs at the end
		print (str(hit[13]) + '\t' + str(hit[12]) + '\t' + hit[1] + '\t' + hit[2] + '\t' + hit[3] + '\t' + str(hit[4]) + '%\t' + str(hit[5]) + '\t' + str(hit[7]) + '\t' + str(hit[14]))

if __name__ == '__main__':
	main()
