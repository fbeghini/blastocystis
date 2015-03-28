#!/usr/bin/env python

import gzip, bz2, os, argparse, sys, sqlite3, csv
from Bio import SeqIO, Entrez
from bz2 import BZ2File

parser = argparse.ArgumentParser()
parser.add_argument("--dataset", "-d", help="Location of the dataset")
parser.add_argument("--output","-o", help='Database files')
parser.add_argument("lsumap")
parser.add_argument("ssumap")

args = parser.parse_args()

Entrez.email	=	"fbeghini (at) live (dot) com"
retmax		=	1000

# Genes: Build index with faidx and associate id -> organism
# Full genomes idem
# Repophlan idem
# rRNA: extract last field from taxonomy
# Genomes: How?
if args.dataset and args.output:

	# try:
	# 	conn = sqlite3.connect(args.database)
	# except:
	# 	print "Unable to connect to "+args.database

	# conn.row_factory = sqlite3.Row
	# cursor = conn.cursor()

	ssu = {}
	lsu = {}

	with open(args.lsumap) as rnamap:
		lsu = dict(csv.reader(rnamap, delimiter='\t')) 
	with open(args.ssumap) as rnamap:
		ssu = dict(csv.reader(rnamap, delimiter='\t')) 

	# cursor.execute("CREATE TABLE IF NOT EXISTS pangenomeMap (id varchar(255), organism VARCHAR(255), PRIMARY KEY(id))")
	if os.path.isdir(args.dataset):
		with open(args.output,"w") as tsv:
			for(path, dirs, files) in os.walk(args.dataset):
				for g in files:
					extention = os.path.splitext(g)[1]
					ids = []
					result = []

					openm = open if any(f in extention for f in ["fasta","fa","fna","ffn"]) else gzip.open if 'gz' in extention else BZ2File if "bz2" in extention else open
					handle = openm(path+"/"+g)
					if 'genes' in path:
						print g
						oname = g.split('/')[-1].split('.')[0]
						for seq in SeqIO.parse(handle, "fasta"):
							tsv.write( "%s;frag_%s;%i\n" % (seq.id, oname, len(seq.seq) ))
					elif 'rRNA' in path:
						print g
						for seq in SeqIO.parse(handle, "fasta"):
							_id = seq.id.split('.')[0]
							type = "28S" if _id in lsu else "18S" if _id in ssu else "rRNA"
							tsv.write( "%s;%s_%s_%s;%i\n" % (seq.id, type, seq.description.split(";")[-1], seq.id, len(seq.seq)))
					elif 'full_genomes' in path:
						print g
						for seq in SeqIO.parse(handle, "fasta"):
							tsv.write( "%s;refg_%s;%i\n" % (seq.id, " ".join(seq.description.split(" ")[1:3]), len(seq.seq)))
					elif 'genomes' in path:
						print g
						for seq in SeqIO.parse(handle, "fasta"):
							ids.append(seq.id.split("|")[1])
						query_len = len(ids)
						print "Len of query is %s" % query_len

						# try:
						for start in range(0, query_len, retmax):
							end = min(query_len,start+retmax)
							query = ','.join(ids[start:end])
							out = open("/tmp/query","w")
							out.write(",".join(ids))
							print "\tFetching names from %i to %i" % (start+1, end)
							handle = Entrez.epost(db="nucleotide", id=query)
							postquery = Entrez.read(handle)
							query_key = postquery["QueryKey"]
							webenv = postquery["WebEnv"]
							handle = Entrez.efetch(db="nucleotide", retmode="xml", query_key=query_key, webenv=webenv)
							result += Entrez.read(handle, validate=False)
							handle.close()
						
						print "Len of result is %s" % len(result)
						for r in result:
							seqids = r["GBSeq_other-seqids"]
							seqids.reverse()
							tsv.write( "%s;genes_%s;%s\n" % ("|".join(seqids), r["GBSeq_source"], r["GBSeq_length"]))
						# except:
						# 	print sys.exc_info()
# NCBI 
