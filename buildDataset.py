#! /usr/bin/env python
import os, argparse, gzip, bz2	
from Bio import SeqIO

parser = argparse.ArgumentParser()
parser.add_argument("input-folder", help="Location of the dataset")
parser.add_argument("name", help="Name of the output file")

args = parser.parse_args()

#Build fasta dataset
if(os.path.isdir(args.input_folder)):
	with open(name,"a") as out:
		for (p,d,f) in os.walk(args.input_folder):
			ext = os.path.splitext(f)[1]
            if any(f in extention for f in ["fasta","fa","fna"]):
		        handle = open(path+"/"+g,"r")
			elif 'gz' in extention:
		        handle = gzip.open(path+"/"+g)
			elif "bz2" in extention:
		        handle = BZ2File(path+"/"+g)

		  	out.write(handle.readlines())

#Convert dataset

SeqIO.convert(name,"fasta",name+".conv","fasta")
os.remove(name)
os.rename(name+".conv",name)
try:
	SeqIO.convert("","fasta","","fasta")
except:
	pass