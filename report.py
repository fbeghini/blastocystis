#! /usr/bin/env python

import csv, argparse, os

parser = argparse.ArgumentParser()

parser.add_argument("TSV")
parser.add_argument("out")

args = parser.parse_args()

report = {}
inputTSV = args.TSV.split(',')

for sample in inputTSV:
	with open(sample) as f:
		reader = csv.reader(f,delimiter='\t')
		reader.next()
		for row in reader:
			organism = row[0]
			coverage = row[1]
			
			if organism not in report:
				report[organism] = dict((k,0) for k in inputTSV)
			report[organism][sample] = coverage


with open("/tmp/out.txt","w") as outfile:
	outfile.write("#%s\n" % "\t".join(inputTSV))
	for organism in report:
			outfile.write("%s\t%s\n" %(organism, '\t'.join([str(x) for x in report[organism].values()])))