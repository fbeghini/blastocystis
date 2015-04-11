#!/usr/bin/env python
import os, csv, argparse, pickle, sys, time, pandas
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
		for row in tsv:
			fields = row.split("\t")		# fieldnames=(id, cov, base, tot, %)
			if fields[0] != "genome":
				try:
					read = fields[0].split("|")[1] if "|" in fields[0] else fields[0]
					contigs[read] += Counter({ int(fields[1]) : int(fields[2]) } )
				except:
					print fields[0] +" is missing, ignoring...Try rebuilding the genome map with buildGenomeMap"

	with open(args.inputTSV.replace(".tsv",".bed")) as bed:
		for row in bed:
			fields = row.split("\t")		# fieldnames=("id", "start", "end", "coverage"))
			if "|" in fields[0]:
				idRead = genomeMap[ next(x for x in genomeMap if fields[0].split("|")[1] in x) ][0]
			else:
				idRead = genomeMap[fields[0]][0]	
			coverageInterval = int(fields[2]) - int(fields[1])
			coverageLevel = float(fields[3])

			oldTot = meta[idRead][0]
			oldCov = meta[idRead][1]
			meta[idRead] = [oldTot+coverageInterval, oldCov+coverageLevel]
	
	for idc, value in contigs.iteritems():
		genome = genomeMap[ next(x for x in genomeMap if idc in x) ][0]
		for depth, coverage in value.iteritems():
			mappedcontigs[genome] += Counter( {depth : coverage} )
	
	for k,v in genomeMap.itervalues():
		glength[k] += v

	with open(args.inputTSV.replace(".tsv","_mappedd.tsv"),"w") as tsv:
		tsv.write('Genome\t%-coverage\tX-coverage\tBase covered\tTotal\n')
		for genome in mappedcontigs:
			tot = glength[genome]
			baseCovered=meta[genome][0]
			totalRead =meta[genome][1]
			pcov = float(baseCovered)/tot
			xcov = float(baseCovered*totalRead)/tot
			tsv.writelines("%s\t%1.4f\t%1.4f\t%i\t%i\n" % (genome, pcov, xcov, baseCovered, tot))