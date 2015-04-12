#!/usr/bin/env python
import argparse
import pandas as pd
import numpy as np

parser = argparse.ArgumentParser()

parser.add_argument("inputBED", help="BED to be mapped")
parser.add_argument("contigGenomeCSV", help="Map contig-genome")

args = parser.parse_args()

if args.inputBED and args.contigGenomeCSV:
	
	genomeMap = pd.read_csv(args.contigGenomeCSV, sep=';', names=['id','org','len'])
	#genomeMap=genomeMap.drop_duplicates('id')
	#genomeMap.to_csv('/tmp/pangenome.csv', sep=';', header=False, index=False)
	genomeMap.id = genomeMap.id.apply(lambda x: x.split("|")[1] if "|" in x else x)
	
	# tsv = pd.read_csv(args.inputTSV, sep=';', names=['id','org','len'])
	# tsv = tsv[tsv['id']!='genome']
	# tsv.id = tsv.id.apply(lambda x: x.split("|")[1] if "|" in x else x)
	
	# contigs = pd.merge(tsv, genomeMap, how='inner', on='id').groupby(['org','cov']).agg({'base':np.sum})

	bed = pd.read_csv(args.inputBED, sep='\t', names=['id', 'start', 'end', 'depthCoverage'])
	bed.id = bed.id.apply(lambda x: x.split("|")[1] if "|" in x else x)
	bed = pd.merge(bed,genomeMap, how='inner', on='id')[['org','start','end','depthCoverage']]
	meta = pd.DataFrame()

	for id, ds in bed.groupby('org'):
		ds['baseCovered'] = ds['end']-ds['start']
		x = ds[['baseCovered','depthCoverage']].sum()
		meta = meta.append({ 'org' : id, 'baseCovered' : x['baseCovered'], 'depthCoverage' : x['depthCoverage'] },ignore_index=True)
	glength = genomeMap.groupby('org', as_index=False).agg({'len':np.sum})

	out=pd.merge(glength, meta)
	out['%']=out.baseCovered/out.len
	out['fold']=(out.depthCoverage*out.baseCovered)/out.len

	out.to_csv(args.inputBED.replace(".bed","_mapped.bed"), sep='\t', index=False, float_format='%1.4f')