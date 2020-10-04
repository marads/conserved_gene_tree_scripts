import sys, os , glob, argparse

def init():
	parser = argparse.ArgumentParser(description='Script to select fastas that have all species represented once and concatenate the alignment')
	parser.add_argument("-inDir", dest='inDir', required=True, help='dirname containing the (aligned and trimmed) fastafiles (end with .fasta) to concatenate')
	parser.add_argument("-outDir", dest='outDir', help='outdirname if output should not be saved in same dir as inDir')
	parser.add_argument("-splitAt", dest='split_str', default = '__', help='Symbol where to split the fasta header, default is "__".')
	parser.add_argument("-nSpecies", dest='number_of_species', required=True, type=int, help='Number of species that should be represented in the multiple sequence alignment')
	
	parser.add_argument("-v", "--verbose", default = False, help="print progress messages", action="store_true")
	

	########################
	#
	#  Catch possible errors here
	#
	########################
	args = parser.parse_args()
	if args.outDir == None:
		args.outDir = args.inDir

	report = '\n\n##################################\n#\n#   SETTINGS\n'
	argsdict = vars(args)
	for var in argsdict.keys():
		report += '#\t'+var+'\t'+str(argsdict[var])+'\n'
	report += '#\n##################################\n\n'
	print(report)

	return args




	

	
def fasta2dicts(fasta_fname, split_str):
	'''
	Reads a fastafile and creates a dictionary: header -> sequence
	Input
	fasta_fname : name of the fasta_file
	split_str   : the string at which to split te header; the character that 
	              separates the species name from the rest of the sequence id.
	Returns
	A dictionary: species name --> sequence
	'''
	id2seq    = {}
	fastafile = open(fasta_fname)
	line      = fastafile.readline()
	assert len(line) > 2, '*** Error! this fastafile is empty!'
	assert line[0] == '>', '*** Error! this file is not in fasta format, I expect a ">" on the first line!'
	
	currentid = line.strip().split(split_str)[0][1:]
	
	seq  = ''
	line = fastafile.readline()
	while len(line)>0:
		if line[0]=='>': # new header
			assert currentid not in id2seq, "Error! key"+currentid+"already exists!"
			
			if seq[-1] == '*': seq = seq[:-1]
			id2seq[currentid]=seq.upper()
		
			# initialize new sequence identifier
			currentid = line.strip().split(split_str)[0][1:]
			
			seq=''
		else:
			# add this line to sequence
			seq += line.strip()
			
		line = fastafile.readline()
	
	#the last one
	if seq[-1] == '*': seq = seq[:-1]
	id2seq[currentid]=seq.upper()
	
	fastafile.close()
	return id2seq



	

def select_and_concatenate(fasta_indir, split_str, number_of_species, outfname):
	'''
	Selects from all fastafiles in <fasta_indir> those that have sequences for 
	<number_of_species> species. Species names are separated from the rest of 
	the sequence identifier in the header by <split_str>.
	Input:
	fasta_indir       : the folder that contains the fasta sequences (ending with '.fasta)
	split_str         : the string at which to split te header; the character that 
	                    separates the species name from the rest of the sequence id.
	number_of_species : the exact number of species in the dataset that should be represented in each fasta file
	outfname          : the name of the file that will contain the concatenated sequences in fasta format.     
	
	Returns:
	The number of fastafiles selected for the concatenated alignment.
	'''

	number_of_concatenated  = 0 
	species2concatenatedseq = {}
	species_set = set([])

	for fasta_fname in glob.glob(fasta_indir+'*.fasta'):

		id2seq = fasta2dicts(fasta_fname, split_str)
		if len(id2seq) == number_of_species:
			species_set_f = set(id2seq.keys())
			if len(species_set) == 0: 
				# initialize set of species
				species_set = species_set.union(species_set_f)
				species2concatenatedseq = id2seq
				number_of_concatenated += 1
				print(fasta_fname)
			elif species_set_f == species_set:
				number_of_concatenated += 1
				lengths = set([])
				for species in id2seq.keys():
					seq = id2seq[species]
					lengths.add(len(seq))
					species2concatenatedseq[species] += seq
				print(fasta_fname)
				assert len(lengths) == 1, 'Error! Not all sequences are the same length while they should be!\n***\t'+fasta_fname+'\n'+str(lengths)
		else:
			print(fasta_fname, len(id2seq))
		

	lengths = set([])
	for species in species2concatenatedseq:
		lengths.add(len(species2concatenatedseq[species]))

	assert len(lengths) == 1, 'Error! Not all sequences are the same length while they should be!'

	with open (outfname, 'w') as outfile:
		for species in species2concatenatedseq:
			seq = species2concatenatedseq[species]
			outfile.write('>'+species+'\n'+seq+'\n')
		
		return number_of_concatenated






if __name__ == "__main__":

	args = init()

	outfname = args.outDir + '/concatenated_alignment2.fasta'
	number_of_concatenated = select_and_concatenate(args.inDir, args.split_str, args.number_of_species, outfname)


			
	# One line, after selecting suitable fasta files
	# RS : record separator
	# NF : number of fields
	# 
	# not necessarily faster
	# cut -f1 -d "__" *.fasta | awk -v RS=">" -v FS="\n" -v OFS="\n" '{if (NF=149) for(i=2; i<=NF; i++) {seq[$1] = seq[$1]$i}}; END {for(id in seq){print ">"id, seq[id]}}' > combined.fa
		

