#!/usr/bin/env python
import os, csv, argparse, pickle, sys, time
from collections import defaultdict, Counter

parser = argparse.ArgumentParser()

parser.add_argument("inputTSV", help="TSV to be mapped")
parser.add_argument("contigGenomeCSV", help="Map contig-genome")

args = parser.parse_args()

if args.inputTSV and args.contigGenomeCSV:
	genomeMap = (lambda: defaultdict(lambda: ('',0)))()

	struct = lambda: defaultdict(Counter)
	_meta = lambda: defaultdict(lambda: [0,0])
	mappedcontigs = struct()
	contigs = struct()                                                              
	meta = _meta()

	glength = Counter()

	with open(args.contigGenomeCSV) as csvfile:
		for item in csv.DictReader(csvfile, delimiter=';', fieldnames=("id","organism","length")):
			genomeMap[item['id']] = (item['organism'], int(item['length']))

	with open(args.inputTSV) as tsv:
		reader = csv.reader(tsv, delimiter='\t')
		for row in reader:
			if row[0] != "genome":
				try:
					read = row[0].split("|")[1] if "|" in row[0] else row[0]
					contigs[read] += Counter({ int(row[1]) : int(row[2]) } )
				except:
					print row[0] +" is missing, ignoring...Try rebuilding the genome map with buildGenomeMap"
	
	# pickle.dump(contigs, open("/home/francesco.beghini/contig.pk","wb"))

	# contigs = pickle.load(open("/home/francesco.beghini/contig.pk","rb"))

	with open(args.inputTSV.replace(".tsv",".bed")) as bed:
		reader = csv.DictReader(bed, delimiter='\t', fieldnames=("id", "start", "end", "coverage"))
		for row in reader:
			if "|" in row["id"]:
				idRead = genomeMap[ next(x for x in genomeMap if row["id"].split("|")[1] in x) ][0]
			else:
				idRead = genomeMap[row["id"]][0]
			coverageInterval = int(row["end"]) - int(row["start"])
			coverageLevel = int(row["coverage"])

			oldTot = meta[idRead][0]
			oldCov = meta[idRead][1]
			meta[idRead] = [oldTot+coverageInterval, oldCov+coverageLevel]
	
	# pickle.dump(meta, open("/home/francesco.beghini/meta.pk","wb"))
	# meta = pickle.load(open("/home/francesco.beghini/meta.pk","rb"))
	
	for idc, value in contigs.iteritems():
		genome = genomeMap[ next(x for x in genomeMap if idc in x) ][0]
		for depth, coverage in value.iteritems():
			mappedcontigs[genome] += Counter( {depth : coverage} )
	
	for k,v in genomeMap.itervalues():
		glength[k] += v

	# pickle.dump(mappedcontigs, open("/home/francesco.beghini/mappedcontigs.pk","wb"))

	# mappedcontigs = pickle.load(open("/home/francesco.beghini/mappedcontigs.pk","rb"))

	with open(args.inputTSV.replace(".tsv","_mapped.tsv"),"w") as tsv:
		tsv.write('Genome\t%-coverage\tX-coverage\tBase covered\tTotal\n')
		for genome in mappedcontigs:
			tot = glength[genome]
			baseCovered=meta[genome][0]
			totalRead =meta[genome][1]
			pcov = float(baseCovered)/tot
			xcov = float(baseCovered*totalRead)/tot
			tsv.writelines("%s\t%1.4f\t%1.4f\t%i\t%i\n" % (genome, pcov, xcov, baseCovered, tot))