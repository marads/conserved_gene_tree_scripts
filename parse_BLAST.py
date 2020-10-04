description = '''
   Here we parse the output from BLAST, where we select homologs --> blast hits that meet our criteria
'''

import argparse
import os
import sys
import glob



def parse_cmndline_arguments():
	'''
	Parse arguments and perform necessary checks on the arguments that were passed (e.g. whether a file actually exists).
	'''
	parser = argparse.ArgumentParser(description = description)
	parser.add_argument('-queryFasta', dest='query_fasta', required=True, type=str, help="Name of the fastafile with the query sequences.")
	parser.add_argument('-dbFasta', dest='db_fasta', required=True, type=str, help="Name of the fastafile with the target/subject/database sequences.")
	parser.add_argument('-BLASToutput', dest='blast_output', required=True, type=str, help="Name of the fastafile with the target/subject/database sequences.")
	parser.add_argument('-Ngenomes', dest='n_genomes', required=True, type=int, help="Number input genomes.")
	parser.add_argument('-outDir', dest='out_dir', required=True, type=str, \
	help="The name of the output directory. Output will be stored in: outDir/02.bedfilePerQuery, outDir/02.fastaPerQuery")
	parser.add_argument('-minOverlap', dest='min_overlap', default = 0.7, type=float, \
	help="Number between 0 and 1 that indicates the fraction of the query that has to be in the HSP")
	parser.add_argument('-minPercIdentity', dest='min_perc_identity', default = 80.0, type=float, \
	help="Number between 0 and 100 that indicates the minimum percent identity between the query and the subject in the HSP")
	parser.add_argument('-maxEvalue', dest='max_evalue', default = 0.001, type=float, help="Maximum E-value")

	#parser.add_argument('--TEST', dest='test', default = False, action = 'store_true', help="Set this flag to test this code using the files in the 'test' folder.")
	args = parser.parse_args()

	#check if files exist
	assert os.path.exists(args.query_fasta), 'The fastafile '+args.query_fasta+' does not exist, please check the filename and the path again.'
	assert os.path.exists(args.db_fasta), 'The fastafile '+args.db_fasta+' does not exist, please check the filename and the path again.'
	assert os.path.exists(args.blast_output), 'The fastafile '+args.blast_output+' does not exist, please check the filename and the path again.'
	assert os.path.exists(args.out_dir), 'The folder '+args.out_dir+' does not exist, please check the path.'

	#create output folder
	os.system('mkdir -p '+args.out_dir+'/02.bedfilePerQuery')
	os.system('mkdir -p '+args.out_dir+'/02.fastaPerQuery')

	report = '\n\n##################################\n#\n#   SETTINGS\n#'
	argsdict = vars(args)
	for var in argsdict.keys():
		report += '#\t'+var+'\t'+str(argsdict[var])+'\n'
	report += '#\n##################################\n\n'
	print(report)

	return args


def fasta2dict(fasta_fname):
	'''
	This reads a fastafile and returns a dictionary with the sequence per ID
	'''
	seqid2seq = {}
	
	seqid = None
	for line in open(fasta_fname).readlines():
		if len(line.strip()) > 0:
			if line[0] == '>': #header
				seqid = line[1:].split()[0]
				seqid2seq[seqid] = ''
			else:
				seqid2seq[seqid] += line.strip()

	#last one:
	seqid2seq[seqid] += line.strip()

	return seqid2seq



def blastout2dict(args, seqid2seq):
	'''
	Read BLAST file, select hits that meet requirements wrt percent identity and overlap etc.
	Puts this data in a nested dictionary: query -> subject -> HSP data
	'''	
	
	query2subject2hits = {}
	
	with open(args.blast_output) as blast_file:
		for line in blast_file.readlines():
			if line[0] != '#': 
				tabs    = line.strip().split('\t')
				query   = tabs[0]
				subject = tabs[1] 
				perc_id = float(tabs[2])
				length  = int(tabs[3])
				qstart, qend, sstart, send = [int(i) for i in tabs[6:-2]]
				evalue  = float(tabs[10])
				qlen    = float(len(seqid2seq[query]))
				
				# test if hit meets requirements: if the hit is close and complete enough to include in a MSA
				if evalue <= args.max_evalue and length/qlen >= args.min_overlap and perc_id >= args.min_perc_identity:
					if not query in query2subject2hits:
						query2subject2hits[query] = {}

					if subject not in query2subject2hits[query]:
						query2subject2hits[query][subject] = []

					query2subject2hits[query][subject].append((perc_id, length, qstart, qend, sstart, send, evalue))

	return query2subject2hits



def blastout2HSPbed(query2subject2hits, Ngenomes, out_dir):
	'''
	Parse the dictionary with BLAST output. Select queries that have a hit (but not more than one) in all subjects and for these queries create a bedfile with hits.
	Input: dictionary query -> subject -> data HSPs
	Output: dictionary query -> name of the bedifle with hits for this query
	'''
	
	query2bed_fname = {}
	for query in query2subject2hits: 
		#print query
		bed_fname = out_dir+'/02.bedfilePerQuery/'+query.split('(')[0]+'.bed'
		
		bedstr = '' # collect info on locations of HSPs in this string, in bed format:
		
		# check if this query has a hit in all subjects 
		if len(query2subject2hits[query]) == Ngenomes:
		# check if this query has only one hit in all subjects: 
			one_hit = True
			for subject in query2subject2hits[query]:
				if len(query2subject2hits[query][subject]) > 1:
					one_hit = False
					break
			
			if one_hit:
				for subject in query2subject2hits[query]:
					for perc_id, length, qstart, qend, sstart, send, evalue in query2subject2hits[query][subject]:
				
						# construct line for the bedfile:
						# The three obligatory fields: sequence id, start, end
						# and three optional fields: name, score, strand
						# Note: sstart-1 because BLAST starts counting at 1, while bedtools starts counting at 0
						bedstr += subject+'\t'
						if sstart < send:
							bedstr += str(sstart-1)+'\t'+str(send)+'\t'+subject+'__'+str(sstart-1)+'-'+str(send)+'\t'+str(evalue)+'\t+\n' 
						else:
							bedstr += str(send-1)+'\t'+str(sstart)+'\t'+subject+'__'+str(sstart)+'-'+str(send-1)+'\t'+str(evalue)+'\t-\n'
		
		
				with open(bed_fname, 'w') as bedfile:
					bedfile.write(bedstr)
				
				query2bed_fname[query] = bed_fname

	return query2bed_fname



def bed2fasta(query2sequence, query2bed_fname, database_fasta):
	'''
	For each query, create a fastafile with the query sequence and the sequence of the HSPs
	database_fasta is the concatenated fasta with the genome sequences
	'''
	for query in query2sequence:
		if query in query2bed_fname: # query has single hit in all genomes
			# create output fastafile
			bed_fname   = query2bed_fname[query]
			fasta_fname = bed_fname.replace('02.bedfilePerQuery','02.fastaPerQuery').replace(".bed", ".fasta")
			
			# write the sequences of the hits to a fasta file
			cmnd = 'bedtools getfasta -fi '+database_fasta+' -bed '+bed_fname+' -s -fo '+fasta_fname
			print(cmnd, os.system(cmnd))
			
			# then add the query sequence
			with open(fasta_fname, 'a') as fasta_file:
				fasta_file.write('>'+query+'\n')
				fasta_file.write(query2sequence[query]+'\n')
				


if __name__ == "__main__":

	# parse arguments passed by the user
	args = parse_cmndline_arguments()
	
	# sequence per query
	seqid2seq = fasta2dict(args.query_fasta)
	
	# blast hits per query, per subject, filtered for overlap etc.
	query2subject2hits = blastout2dict(args, seqid2seq)
	
	# bedfile per query, queries filtered for having a single hit in all genomes
	query2bedfile = blastout2HSPbed(query2subject2hits, args.n_genomes, args.out_dir)
	
	# fasta file per query
	bed2fasta(seqid2seq, query2bedfile, args.db_fasta)
	
	