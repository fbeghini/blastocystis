#! /usr/bin/env python
import os, argparse, gzip, removeDuplicates
from Bio import SeqIO
from bz2 import BZ2File
from tempfile import NamedTemporaryFile

#Build fasta dataset
def buildDataset(inputFolder, outName):
	with NamedTemporaryFile(delete=True) as _out:
		for (path,dirs,files) in os.walk(inputFolder):
			for f in files:
				extension = os.path.splitext(f)[1]
				if any(f in extension for f in ["fasta","fa","fna"]):
					handle = open
				elif 'gz' in extension:
					handle = gzip.open
				elif "bz2" in extension:
					handle = BZ2File
				with handle(path+"/"+f) as genome:
					_out.writelines(genome.readlines())
		SeqIO.convert(_out.name,"fasta",outName,"fasta")

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument("inputfolder", help="Location of the dataset")
	parser.add_argument("name", help="Name of the output file")

	args = parser.parse_args()
	if(os.path.isdir(args.inputfolder)):
		buildDataset(args.inputfolder, args.name)
		removeDuplicates.filter(args.name)