#!/usr/bin/env python

import gzip, bz2, os, argparse, sys, sqlite3
from Bio import SeqIO, Entrez
from bz2 import BZ2File

parser = argparse.ArgumentParser()
parser.add_argument("--dataset", "-d", help="Location of the dataset")
parser.add_argument("--database","-o", help='Database files')

args = parser.parse_args()

Entrez.email	=	"francesco.beghini@studenti.unitn.it"
retmax			=	10000
db				=	"nucleotide"
retmode 		=	"xml"

# Genes: Build index with faidx and associate id -> organism
# Full genomes idem
# Repophlan idem
# rRNA: extract last field from taxonomy
# Genomes: How?
if args.dataset and args.database:

	# try:
	# 	conn = sqlite3.connect(args.database)
	# except:
	# 	print "Unable to connect to "+args.database

	# conn.row_factory = sqlite3.Row
	# cursor = conn.cursor()

	# cursor.execute("CREATE TABLE IF NOT EXISTS pangenomeMap (id varchar(255), organism VARCHAR(255), PRIMARY KEY(id))")
	if os.path.isdir(args.dataset):
		for(path, dirs, files) in os.walk(args.dataset):
			with open(args.output,"a+") as tsv:
				for g in files:
					extention = os.path.splitext(g)[1]
					ids = []
					result = []

					if any(f in extention for f in ["fasta","fa","fna"]):
						handle = open(path+"/"+g,"r")
					elif 'gz' in extention:
						handle = gzip.open(path+"/"+g)
					elif "bz2" in extention:
						handle = BZ2File(path+"/"+g)

					if 'genes' in path:
						print path, g
						oname = g.split('/')[-1].split('.')[0]
						for seq in SeqIO.parse(handle, "fasta"):
							tsv.write( "%s;%s\n" % (seq.id, oname))
					elif 'rRNA' in path:
						print path, g0
						for seq in SeqIO.parse(handle, "fasta"):
							tsv.write( "%s;%s\n" % (seq.id, seq.description.split(";")[-1]))
					elif 'full_genomes' in path:
						print path, g
						for seq in SeqIO.parse(handle, "fasta"):
							tsv.write( "%s;%s\n" % (seq.id, " ".join(seq.description.split(" ")[1:3])))
					elif 'genomes' in path:
						print path, g
						for seq in SeqIO.parse(handle, "fasta"):
							ids.append(seq.id.split("|")[1])
						query_len = len(ids)
						query = ','.join(ids)
						print "len of query is %s" % query_len
						try:
							postquery = Entrez.read(Entrez.epost(db="nucleotide", id=query))
							query_key = postquery["QueryKey"]
							webenv = postquery["WebEnv"]
							retmax = 10000

							for retstart in range(0, query_len, retmax):
								handle = Entrez.efetch(db=db, retmode=retmode, retmax=retmax, retstart=retstart, query_key=query_key, webenv=webenv)
								_results = Entrez.read(handle, validate=False)
								result += _results

							print "Len of result is %s" % len(result)
							for r in result:
								seqids = r["GBSeq_other-seqids"]
								seqids.reverse()
								tsv.write( "%s;%s\n" % ("|".join(seqids), r["GBSeq_source"]))
						except:
							print sys.exc_info()
