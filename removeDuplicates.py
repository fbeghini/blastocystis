#!/usr/bin/env python
import sys, argparse
from Bio import SeqIO
from Bio.SeqUtils.CheckSum import seguid

def sequence_cleaner(fastafile):
	ids = set()
	checksums=set()
	for record in SeqIO.parse(fastafile, "fasta"):
		if seguid(record.seq) in checksums or record.id in ids:
			print "Ignoring %s" % record.description
			continue
		checksums.add(seguid(record.seq))
		ids.add(record.id)
		yield record

def filter(inputFASTA):
	records = sequence_cleaner(inputFASTA)
	count = SeqIO.write(records, inputFASTA.replace(".fasta",".cleared.fasta"), "fasta")
	return count

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument("inputFASTA", metavar="inputFASTA", help="FASTA file to be cleared")
	args = parser.parse_args()

	try:
		count = filter(args.inputFASTA)
		print "Saved %d records" % count
	except:
		print sys.exc_info()[0]