from Bio import SeqIO
import os, sys, argparse

parser = argparse.ArgumentParser()

groupin = parser.add_mutually_exclusive_group()
groupin.add_argument("--folder", help="Folder with fastas")
groupin.add_argument("--file", help="Single fasta")

parser.add_argument("--pattern", help="Pattern contained in description")
parser.add_argument("-v", help="Filter behaviour, [R]emove the sequences containing the pattern, [S]ave the sequence containing the pattern", choices=["R","S"], default="R")

args = parser.parse_args()

fl=[]
pattern = args.pattern.split(',')

remove =  True if args.v =='R' else False

if args.folder and os.path.isdir(args.folder):														#Folder with multiple fastas
	for (path, dirs, files) in os.walk(args.folder):
		if len(files)>0:
			fl.append([path+"/"+f for f in files][0])
elif args.file and os.path.isfile(args.file):														#Single fasta file
	fl.append(args.file)

try:
	for fastas in fl:
		with open(fastas+".filtered", "w+") as out:
			for record in SeqIO.parse(fastas,"fasta"):
				if(remove):
					if(not any(p.lower() in record.description.lower() for p in pattern)):						#remove
						print "%s:>%s" % (fastas, record.description)
						out.write(">%s\n%s\n" % (record.description, record.seq))
				else:
					if(any(p.lower() in record.description.lower() for p in pattern)):							#save
						print "%s:>%s" % (fastas, record.description)
						out.write(">%s\n%s\n" % (record.description, record.seq))				
except:
	print sys.exc_info()[0]