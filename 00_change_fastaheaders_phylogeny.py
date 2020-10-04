# created 13-07-2019 by Mara de Sain
# changes fastaheaders to >filename__contignumber
# Use: python change_fastaheaders.py [absolute path to directory containing genome .fasta files] 


import os, sys, errno
from Bio import SeqIO

def main(directory):
    outputpath = directory + "/output"
    create_outfolder(outputpath)
    loop_over_fasta_files(directory, outputpath)
    
# make outfolder if it does not exist yet in python 2.7 --> nicer way in python 3: os.makedir(os.path.dirname(outputfolder), exist_ok=True)
def create_outfolder(outputpath):
    print "Outputfolder is: " + outputpath
    if not os.path.exists(outputpath):
        try:
            os.mkdir(outputpath)
        except OSError as exc: # Guard against race condition
            if exc.errno != errno.EEXIST:
                raise

def loop_over_fasta_files(directory, outputpath):
#loop over files in folder
    for filename in os.listdir(directory):
        if filename.endswith(".fasta"):
            name = filename[:-6]
            filepath = directory + "/" + filename
            parse_fasta_files(filename, filepath, name, outputpath)
        else:
            pass

def parse_fasta_files(filename, filepath, name, outputpath):
    with open (filepath, "rU") as fastafile:
        for record in SeqIO.parse(fastafile, "fasta"):
            new_record_id = create_new_record_id(record, name)
            write_new_records_to_file(filename, new_record_id, record, outputpath)                 

def create_new_record_id(record, name):
    split_id = record.id.split("_")
    new_record_id = name + "__contig" + split_id[-1]
    return new_record_id
    
def write_new_records_to_file(filename, new_record_id, record, outputpath):
    outfilepath = outputpath + "/" + filename
    with open (outfilepath, "a+") as outfile:
        new_fasta_sequence = str(">" + new_record_id + "\n" + record.seq + "\n")
        outfile.write(new_fasta_sequence)

main(sys.argv[1])