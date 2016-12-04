#! /usr/bin/env python
import csv, argparse, glob, re, os

parser = argparse.ArgumentParser()

parser.add_argument("BED")

args = parser.parse_args()

preport = {}
xreport = {}
rreport = {}
inputBED = glob.glob(args.BED)
#print inputBED
for sample in inputBED:
	with open(sample) as f:
		reader = csv.reader(f,delimiter='\t')
		reader.next()
		for row in reader:
			organism = row[0]
			xcov = row[4]
			pcov = row[5]
			#relab = row[6]

			if organism not in preport and organism not in xreport:# and organism not in rreport:
				preport[organism] = dict((os.path.split(k)[-1],0) for k in inputBED)
				xreport[organism] = dict((os.path.split(k)[-1],0) for k in inputBED)
				rreport[organism] = dict((os.path.split(k)[-1],0) for k in inputBED)
			preport[organism][os.path.split(sample)[-1]] = pcov
			xreport[organism][os.path.split(sample)[-1]] = xcov
			#rreport[organism][os.path.split(sample)[-1]] = relab


with open("%s/breadth_merged.txt" % (os.path.split(inputBED[0])[0]), "w") as outfile:
	outfile.write("organism\t%s\n" % "\t".join(preport.values()[0].keys()))
	for organism in sorted(preport):
			outfile.write("%s\t%s\n" %(organism, '\t'.join([str(preport[organism][x]) for x in preport[organism]])))

with open("%s/depth_merged.txt" % (os.path.split(inputBED[0])[0]), "w") as outfile:
	outfile.write("organism\t%s\n" % "\t".join(xreport.values()[0].keys()))
	for organism in sorted(xreport):
			outfile.write("%s\t%s\n" %(organism, '\t'.join([str(xreport[organism][x]) for x in xreport[organism]])))

# with open("%s/abundances_merged.txt" % (os.path.split(inputBED[0])[0]), "w") as outfile:
# 	outfile.write("organism\t%s\n" % "\t".join(rreport.values()[0].keys()))
# 	for organism in rreport:
# 			outfile.write("%s\t%s\n" %(organism, '\t'.join([str(rreport[organism][x]) for x in rreport[organism]])))
