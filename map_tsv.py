#!/usr/bin/env python

import os, csv, argparse, pickle
from collections import defaultdict

parser = argparse.ArgumentParser()

parser.add_argument("inputTSV", help="TSV to be mapped")
parser.add_argument("outputTSV", help="Output TSV")
parser.add_argument("contigGenomeCSV", help="Map contig-genome")

args = parser.parse_args()

if args.inputTSV and args.outputTSV and args.contigGenomeCSV:
	map = {}																		# Dict of {read:genome}

	mappedcontigs = defaultdict(dict)
	contigs = defaultdict(dict)														# Dict of dict of list of read of a genome
																					# k:v -> {contig as string : 
																					#			{coverage as int : [contig as int, total as int] }
																					#		 }
	with open(args.contigGenomeCSV,"r") as csvfile:
		map = dict(csv.reader(csvfile, delimiter=';')) 

	with open(args.inputTSV,"r") as tsv:
		reader = csv.reader(tsv, delimiter='\t')
		for row in reader:
			if row[0] != "genome":
				try:
					if "|" in row[0]:	
						read = row[0].split("|")[1]
						genome = map[[x for x in map.keys() if read in x][0]]								# Foreach contig, map contig on genome	
					else:
						read = row[0]
						genome = map[read]
					contigs[read][int(row[1])] = [int(row[2]), int(row[3])]							# new contig added, saved covered and total
				except:
					print row[0]

	# pickle.dump(contigs, open("/home/francesco.beghini/contig.pk","wb"))

	# contigs = pickle.load(open("/home/francesco.beghini/contig.pk","rb"))

	#print contigs.keys()
	for curr in contigs:
		completeid = [x for x in map.keys() if curr in x][0]
		genome = map[completeid]
		for covered in contigs[curr]:	
			if covered in mappedcontigs[genome]:
				oldcov = mappedcontigs[genome][covered][0] 
				oldtot = mappedcontigs[genome][covered][1]
				mappedcontigs[genome][covered] = [oldcov + contigs[curr][covered][0], oldtot + contigs[curr][covered][1]]
			else:
				mappedcontigs[genome][covered] = [contigs[curr][covered][0], contigs[curr][covered][1]]
	# scorro ogni contig, assegno un genoma e sommo i valori

	try:
		with open(args.outputTSV,"w") as tsv:
			for genome in mappedcontigs:
				for contig in mappedcontigs[genome]:
					cov=mappedcontigs[genome][contig][0]
					tot=mappedcontigs[genome][contig][1]
					freq=float(cov)/tot
					tsv.write('%s\t%i\t%i\t%i\t%f\n' % (genome, contig, cov, tot, freq))
	except:
		print sys.exc_info()[0]