#! /usr/bin/env python
import csv, argparse, glob, re

parser = argparse.ArgumentParser()

parser.add_argument("BED")

args = parser.parse_args()

preport = {}
xreport = {}
inputBED = glob.glob(args.BED)
for sample in inputBED:
	with open(sample) as f:
		reader = csv.reader(f,delimiter='\t')
		reader.next()
		for row in reader:
			organism = row[0]
			pcov = row[4]
			xcov = row[5]

			if organism not in preport and organism not in xreport:
				preport[organism] = dict((k,0) for k in inputBED)
				xreport[organism] = dict((k,0) for k in inputBED)
			preport[organism][sample] = pcov
			xreport[organism][sample] = xcov

with open("perc_merged.txt","w") as outfile:
	outfile.write("organism\t%s\n" % "\t".join(preport.values()[0].keys()))
	for organism in preport:
			outfile.write("%s\t%s\n" %(organism, '\t'.join([str(preport[organism][x]) for x in preport[organism]])))

with open("fold_merged.txt","w") as outfile:
	outfile.write("organism\t%s\n" % "\t".join(xreport.values()[0].keys()))
	for organism in xreport:
			outfile.write("%s\t%s\n" %(organism, '\t'.join([str(xreport[organism][x]) for x in xreport[organism]])))