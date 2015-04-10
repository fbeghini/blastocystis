#! /usr/bin/env python
import csv, argparse, os, re

parser = argparse.ArgumentParser()

parser.add_argument("TSV")
parser.add_argument("out")

args = parser.parse_args()

report = {}
inputTSV = re.split(',|\n| |',args.TSV)
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

with open(args.out,"w") as outfile:
	outfile.write("#\t%s\n" % "\t".join(inputTSV))
	for organism in sorted(report):
			outfile.write("%s\t%s\n" %(organism, '\t'.join([str(report[organism][x]) for x in sorted(report[organism])])))