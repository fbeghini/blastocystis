from Bio import SeqIO
import os, sys, argparse

path = sys.argv[1]
pattern = sys.argv[2].split(",")

parser = argparse.ArgumentParser()

group = parser.add_mutually_exclusive_group()
group.add_argument("--folder", help="Folder with fastas")
group.add_argument("--file", help="Single fasta")

parser.add_argument("--pattern", help="Pattern contained in description of sequences to be removed")

args = parser.parse_args()

fl=[]

if args.folder and os.path.isdir(args.folder):														#Folder with multiple fastas
	for (path, dirs, files) in os.walk(args.folder):
		if len(files)>0:
			fl.append([path+"/"+f for f in files][0])
elif args.file and os.path.isfile(args.file):														#Single fasta file
	fl.append(args.file)

try:

	for fastas in fl:
		with open(fastas.replace(".fasta",".fa"), "w+") as out:
			with open(fastas, "r+") as handle:
				for record in SeqIO.parse(handle,"fasta"):
					if(not any(p in record.description for p in pattern)):
						print "%s:>%s" % (fastas, record.description)
						out.write(">%s\n%s\n" % (record.description, record.seq))
except:
	print sys.exc_info()[0]