#!/usr/bin/env python
from Bio import SeqIO
import os, gzip, argparse

ids = {}
out = open("/CIBIO/sharedCM/projects/mycobacteria/silva/rRNAseq.fasta","w")
#gunzip fastas
if(not os.path.isfile("/tmp/SILVA_119_SSURef_Nr99_tax_silva.fasta") or not os.path.isfile("/tmp/SILVA_119_LSURef_tax_silva.fasta") or not os.path.isfile("/tmp/SILVA_119_SSUParc_tax_silva.fasta")):
	os.system("zcat /CIBIO/sharedCM/projects/mycobacteria/silva/SILVA_119_SSURef_Nr99_tax_silva.fasta.gz > /tmp/SILVA_119_SSURef_Nr99_tax_silva.fasta")
	os.system("zcat /CIBIO/sharedCM/projects/mycobacteria/silva/SILVA_119_LSURef_tax_silva.fasta.gz > /tmp/SILVA_119_LSURef_tax_silva.fasta")
	os.system("zcat /CIBIO/sharedCM/projects/mycobacteria/silva/SILVA_119_SSUParc_tax_silva.fasta.gz > /tmp/SILVA_119_SSUParc_tax_silva.fasta")

ssu = SeqIO.index("/tmp/SILVA_119_SSURef_Nr99_tax_silva.fasta", "fasta")
lsu = SeqIO.index("/tmp/SILVA_119_LSURef_tax_silva.fasta", "fasta")
ssp = SeqIO.index("/tmp/SILVA_119_SSUParc_tax_silva.fasta","fasta")

for item in open("/home/francesco.beghini/mycobacteria/rRNA","r").readlines():
	key = 'ssu' if 'ssu' in item.split(":")[0] else 'lsu'	
	if key not in ids.keys():
		ids[key] = []
	ids[key].append(item.split(":")[1].strip())

s_keys = [s for s in ssu.keys() for id in ids['ssu'] if id.strip() in s]
l_keys = [l for l in lsu.keys() for id in ids['lsu'] if id.strip() in l]
p_keys = [p for p in ssp.keys() for id in ids['ssu'] if id.strip() in p]

for id in s_keys:
	out.write(">%s\n%s\n" % (ssu[id].description, ssu[id].seq))
for id in l_keys:
	out.write(">%s\n%s\n" % (lsu[id].description, lsu[id].seq))
for id in p_keys:
	out.write(">%s\n%s\n" % (ssp[id].description, ssp[id].seq))

out.close()

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument("inputBED", help="BED to be mapped")
	parser.add_argument("contigGenomeCSV", help="Map contig-genome")
	args = parser.parse_args()