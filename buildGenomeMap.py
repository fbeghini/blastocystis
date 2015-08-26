#!/usr/bin/env python

import gzip, os, argparse, sys, sqlite3, csv
from Bio import SeqIO, Entrez
from bz2 import BZ2File

# Genes: Build index with faidx and associate id -> organism
# Full genomes idem
# Repophlan idem
# rRNA: extract last field from taxonomy
# Genomes: How?

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument("--dataset", "-d", help="Location of the dataset")
	parser.add_argument("--output","-o", help='Database files')
	parser.add_argument("--silva_lsu")
	parser.add_argument("--silva_ssu")
	parser.add_argument("--repophlan_map")

	args = parser.parse_args()

	Entrez.email	=	"matanula@gmail.com"	#Let's trick NCBI changing email each time
	retmax		=	1000
	if args.dataset and args.output:
		lsu = {}
		ssu = {}
		repophlan = {}

		if args.silva_lsu and args.silva_ssu:
			with open(args.silva_lsu) as rnamap:
				lsu = dict(csv.reader(rnamap, delimiter='\t')) 
			with open(args.silva_lsu) as rnamap:
				ssu = dict(csv.reader(rnamap, delimiter='\t')) 
		if args.repophlan_map:
			with open(args.repophlan_map) as repophlanmap:
				for row in repophlanmap:
					id = row.split('\t')[0]
					org = row.split('\t')[18]
					repophlan[id] = org
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
								tsv.write( "%s;genes_%s_%s;%i\n" % (seq.id, oname, seq.id, len(seq.seq) ))
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
						elif 'repophlan' in path:
							print g
							accession = g.split('/')[-1].split('.')[0]
							org = repophlan[accession]
							for seq in SeqIO.parse(handle, "fasta"):
								tsv.write( "%s;refg_%s;%s\n" % (seq.id, org, len(seq.seq)))
						elif 'genus_genes' in path:
							print g
							for seq in SeqIO.parse(handle, "fasta"):
								try:
									ids.append(seq.id.split("|")[1])
								except:
									continue
							ids = list(set(ids))
							query_len = len(ids)
							print "Len of query is %s" % query_len

							try:
								for start in range(0, query_len, retmax):
									end = min(query_len,start+retmax)
									query = ','.join(ids[start:end])
									print "\tFetching names from %i to %i" % (start+1, end)
									handle = Entrez.epost(db="nucleotide", id=query)
									postquery = Entrez.read(handle)
									query_key = postquery["QueryKey"]
									webenv = postquery["WebEnv"]
									handle = Entrez.efetch(db="nucleotide", retmode="xml", query_key=query_key, webenv=webenv)
									result.extend(Entrez.read(handle, validate=False))
									handle.close()
								
								print "Len of result is %s" % len(result)
								for r in result:
									seqids = r["GBSeq_other-seqids"]
									seqids.reverse()
									tsv.write( "%s;genes_%s_%s;%s\n" % ("|".join(seqids), r["GBSeq_source"], "|".join(seqids), r["GBSeq_length"]))
							except:
							 	print sys.exc_info()