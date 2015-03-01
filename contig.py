import gzip, os, argparse
from Bio import SeqIO
from sys import argv


parser = argparse.ArgumentParser("Extract all ")
parser.add_argument("", metavar="inputTSV", help="TSV to be mapped")

if len(argv)<2:
	print "usage"
else:
	sampleFolder = arv[1]
	
	genomes = os.listdir("/banche_dati/sharedCM/Dropbox (CIBIO)/mycobacteria/data/genomes/samples")
	out = open("contig_sample.csv","w+")
#genomes.append(os.listdir("/banche_dati/sharedCM/Dropbox (CIBIO)/mycobacteria/data/genomes/samples"))

	samples = open("mycobacteria_genomes.csv","r")
	dic = {}
	for l in samples:
		_l = l.split(",")
		dic[_l[0]] = _l[1]

	for g in genomes:
		handle = gzip.open("/banche_dati/sharedCM/Dropbox (CIBIO)/mycobacteria/data/genomes/samples/"+g)
		for seq_record in SeqIO.parse(handle,"fasta"):
			if(dic.get(g.split("-")[0],None)!=None):
				out.write(seq_record.id + ";" + str(dic.get(g.split("-")[0],None))+"\n")