#!/usr/bin/env python
import argparse, os, subprocess
import pandas as pd
import numpy as np

f = lambda x: x.split("|")[1] if "|" in x else x

def reads_sample(bamfile):
	idx = subprocess.Popen("samtools idxstats "+bamfile, shell=True, stdout=subprocess.PIPE)
	bam = {}
	for line in idx.stdout:
		id, len, mr, nr = line.strip().split()
		if(int(mr)!=0):
			bam[f(id)] = int(mr)
	return bam

def calculateRelativeAbb(inputBED, contigGenome):
	totalReadsMetagenome = pd.read_table("/scratch/sharedCM/metagenomes//samplereads", sep='\t', names=['dataset', 'sample_name', 'reads'])

	genomeMap = pd.read_csv(contigGenome, sep=';', names=['id','org','len'])
	genomeMap.id = genomeMap.id.apply(lambda x: x.split("|")[1] if "|" in x else x)
	genomeMap=genomeMap.drop_duplicates('id')

	mappedBED = pd.read_table(inputBED)
	mappedBED = mappedBED[mappedBED.org.str.match("ST.$|ST4_WR1")]
	bamReads = reads_sample(inputBED.replace("_mapped.bed",".bam"))
	bam = pd.DataFrame()
	bam['id'] = bamReads.keys()
	bam['reads'] = bamReads.values()
	bam = pd.merge(bam,genomeMap)
	bam = bam[bam.org.str.match("ST.$|ST4_WR1")]
	bam = bam.groupby('org').agg({'reads': np.sum})
	mappedBED.loc[:,"totalReads"] = bam.reads.tolist()
	try:
		totalReads=totalReadsMetagenome.query('sample_name=="'+ inputBED.split('_mapped.bed')[0] + '"').reads.item()#'" & dataset=="Chatelier_gut_obesity"').reads.item()
		mappedBED["relativeAbundance"] = mappedBED.totalReads / totalReads
		mappedBED.to_csv(inputBED, sep='\t', index=False)
	except:
		print("Cannot map " + inputBED)

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument("inputBED", help="BED to be mapped")
	parser.add_argument("contigGenomeCSV", help="Map contig-genome")
	args = parser.parse_args()
	if os.stat(args.inputBED).st_size > 0:
		calculateRelativeAbb(args.inputBED, args.contigGenomeCSV)
 
