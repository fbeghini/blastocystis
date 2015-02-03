import gzip, os
from Bio import SeqIO


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
