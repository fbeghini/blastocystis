#!/usr/bin/env python
import os, csv, argparse, pickle, sys
from collections import defaultdict

parser = argparse.ArgumentParser()

parser.add_argument("inputTSV", help="TSV to be mapped")
parser.add_argument("contigGenomeCSV", help="Map contig-genome")

args = parser.parse_args()

if args.inputTSV and args.contigGenomeCSV:
	genomeMap = {}																		# Dict of {read:genome}
	meta = {}	


	mappedcontigs = defaultdict(dict)
	contigs = defaultdict(dict)														# Dict of dict of list of read of a genome
																					# k:v -> {id as string : 
																					#			{basecovered as int : contig as int }
																					#		 }
	with open(args.contigGenomeCSV) as csvfile:
		genomeMap = dict(csv.reader(csvfile, delimiter=';')) 

	with open(args.inputTSV) as tsv:
		reader = csv.reader(tsv, delimiter='\t')
		for row in reader:
			if row[0] != "genome":
				try:
					read = row[0].split("|")[1] if("|" in row[0]) else row[0]
					contigs[read][int(row[1])] = int(row[2])								# new contig added, saved depth and base covered
				except:
					print row[0] +" is missing, ignoring...Try rebuilding the genome map with buildGenomeMap"

	#pickle.dump(contigs, open("/home/francesco.beghini/contig.pk","wb"))

	# contigs = pickle.load(open("/home/francesco.beghini/contig.pk","rb"))
	
	with open(args.inputTSV.replace(".tsv",".bed")) as bed:
		reader = csv.reader(bed, delimiter='\t')
		for row in reader:
			if "|" in row[0]:
				idRead = genomeMap[ next(x for x in genomeMap if row[0].split("|")[1] in x) ]
			else:
				idRead = genomeMap[row[0]]
			coverageInterval = int(row[2]) - int(row[1])
			coverageLevel = int(row[3])

			if idRead not in meta:
				meta[idRead] = [0,0]
			oldTot = meta[idRead][0]
			oldCov = meta[idRead][1]
			meta[idRead] = [oldTot+coverageInterval, oldCov+coverageLevel]


	for idc in contigs:
		genome = genomeMap[ next(x for x in genomeMap if idc in x) ]
		for covered in contigs[idc]:
			if covered in mappedcontigs[genome]:
				oldpcov = mappedcontigs[genome][covered]
				mappedcontigs[genome][covered] = oldpcov + contigs[idc][covered]			# update base coverage
			else:
				mappedcontigs[genome][covered] = contigs[idc][covered]						# save base covered
				# if idc == "239504781":
				# 	print mappedcontigs[genome]
	# scorro ogni contig, assegno un genoma e sommo i valori

	# pickle.dump(mappedcontigs, open("/home/francesco.beghini/mappedcontigs.pk","wb"))

	# mappedcontigs = pickle.load(open("/home/francesco.beghini/mappedcontigs.pk","rb"))
	

	# try:
	with open(args.inputTSV.replace(".tsv","_mapped.tsv"),"w") as tsv:
		tsv.write('Genome\t%-coverage\tX-coverage\tBase covered\tTotal\n')
		for genome in mappedcontigs:
			tot = sum([mappedcontigs[genome][x] for x in mappedcontigs[genome]])
			baseCovered=meta[genome][0]
			totalRead =meta[genome][1]
			pcov = float(baseCovered)/tot
			xcov = float(baseCovered*totalRead)/tot
			tsv.writelines("%s\t%1.4f\t%1.4f\t%i\t%i\n" % (genome, pcov, xcov, baseCovered, tot))
	# except:
		# print sys.exc_info()[0]