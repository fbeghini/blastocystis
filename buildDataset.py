#! /usr/bin/env python
import os, argparse, gzip, removeDuplicates, subprocess
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

def createBTindex(dataset):
	os.mkdir(dataset)
	build = subprocess.Popen("bowtie2-build %s %s/%s" % (dataset, dataset, dataset), shell=True)
	build.wait()

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument("inputFolder", help="Location of the dataset")
	parser.add_argument("outputFASTA", help="Name of the output file")

	args = parser.parse_args()
	if(os.path.isdir(args.inputfolder)):
		buildDataset(args.inputfolder, args.outputFASTA)
		removeDuplicates.filter(args.outputFASTA)
		createBTindex(args.outputFASTA.replace(".fasta",".cleared.fasta"))
		print "Bowtie2 index built. Location: ./%s" % (args.outputFASTA)