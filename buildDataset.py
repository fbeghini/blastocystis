#! /usr/bin/env python
import os, argparse, gzip  
from Bio import SeqIO
from bz2 import BZ2File
from tempfile import NamedTemporaryFile

parser = argparse.ArgumentParser()
parser.add_argument("inputfolder", help="Location of the dataset")
parser.add_argument("name", help="Name of the output file")

args = parser.parse_args()

#Build fasta dataset
if(os.path.isdir(args.inputfolder)):
	with NamedTemporaryFile(delete=True) as _out:
		for (path,dirs,files) in os.walk(args.inputfolder):
			for f in files:
				extention = os.path.splitext(f)[1]
				if any(f in extention for f in ["fasta","fa","fna"]):
						handle = open(path+"/"+f,"r")
				elif 'gz' in extention:
						handle = gzip.open(path+"/"+f)
				elif "bz2" in extention:
						handle = BZ2File(path+"/"+f)
				_out.writelines(handle.readlines())

		SeqIO.convert(_out.name,"fasta",args.name,"fasta")