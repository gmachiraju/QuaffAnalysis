#! usr/bin/python

#-------------------------------------------------------------------------------------
# This Python script returns splits of data files (FASTQ and FASTA)
# @author: Gautam Machiraju
# Purpose: Holmes Lab 
# Version: Python 2.7.10
#
# Usage: python file_splitter.py <fileName> <number_of_splits>
#-------------------------------------------------------------------------------------
import sys
from itertools import groupby
import textwrap
import re
from Bio import SeqIO

#---------------------------
# Error Checking & Handling
#---------------------------
#Initial checking for a valid file input by trying to open it.
try:
	fileName = sys.argv[1]
	periodLocation = fileName.rfind('.')
	extension = str(fileName)[periodLocation:]
	n = int(sys.argv[2])
	inFile = open(fileName)  

# Exception: if fileName is not able to be opened.
except IOError:
	print "Error: the input file's format could not be opened."
	sys.exit()

# Exception: if user didn't enter a fileName.
except IndexError: 
	print "Error: you must specify: arg1 = fasta/fastq file, arg2 = number of splits." 
	print "Command line should read: 'python file_splitter.py <fileName> <number_of_splits>'"
	sys.exit()

# Exception: if fileName is not located.
except FileNotFoundError: 
	print "Error: the input file you entered could not be opened. Check to make sure it "
	print "exists in the current directory and that you spelled it correctly in the command line."
	sys.exit()

# Checking if file input is of FASTA format (checking for popular extensions).
if (extension == ".fasta") or (extension == ".fa") or (extension == ".fas") or (extension == ".fna") or (extension == ".ffn") or (extension == ".faa") or (extension == ".frn"):
	format = 'fasta'
# Checking if file input is of .txt format.
elif (extension == ".fastq"): 
	format = "fastq"
else:
	print "Error: cannot read file -- it is not of FASTQ or FASTA format."
	sys.exit()


#--------------
# Main Routine
#--------------

# Used for writing files
name = str(fileName)[:periodLocation]

# Collecting read data within file.
records = list(SeqIO.parse(fileName, format))
readNumber = len(records)
totalLength = 0
for record in records:
	totalLength += len(record)

if n > readNumber:
    print "Error: <number_of_splits> is too large."
    sys.exit()

splitLength = readNumber/n


#Adapted from BioPython Website
def batch_iterator(iterator, batch_size) :
    """Returns lists of length batch_size.
 
    This can be used on any iterator, for example to batch up
    SeqRecord objects from Bio.SeqIO.parse(...), or to batch
    Alignment objects from Bio.AlignIO.parse(...), or simply
    lines from a file handle.
 
    This is a generator function, and it returns lists of the
    entries from the supplied iterator.  Each list will have
    batch_size entries, although the final list may be shorter.
    """
    entry = True #Make sure we loop once
    while entry :
        batch = []
        while len(batch) < batch_size :
            try :
                entry = iterator.next()
            except StopIteration :
                entry = None
            if entry is None :
                #End of file
                break
            batch.append(entry)
        if batch :
            yield batch

record_iter = SeqIO.parse(open(fileName),format)
for i, batch in enumerate(batch_iterator(record_iter, splitLength)):

    # Deals with getting 1% of enties
    if n >= 50 and i == 1:
        print "Out of " + str(readNumber) + " total reads"
        sys.exit()

    filename = (name + "_%i" + "_group." + format) % (i+1)
    handle = open(filename, "w")
    count = SeqIO.write(batch, handle, format)
    handle.close()
    print "Wrote %i records to %s" % (count, filename)

print "Out of " + str(readNumber) + " total reads"
inFile.close()




