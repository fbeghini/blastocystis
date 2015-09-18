#!/usr/bin/env python
import argparse, os
import pandas as pd
import numpy as np
import relativeAbundances

def bedmap(inputBED, contigGenome):
	genomeMap = pd.read_csv(contigGenome, sep=';', names=['id','org','len'])
	genomeMap.id = genomeMap.id.apply(lambda x: x.split("|")[1] if "|" in x else x)
	genomeMap=genomeMap.drop_duplicates('id')
	#genomeMap.to_csv('/tmp/pangenome.csv', sep=';', header=False, index=False)
	
	# tsv = pd.read_csv(args.inputTSV, sep=';', names=['id','org','len'])
	# tsv = tsv[tsv['id']!='genome']
	# tsv.id = tsv.id.apply(lambda x: x.split("|")[1] if "|" in x else x)
	
	# contigs = pd.merge(tsv, genomeMap, how='inner', on='id').groupby(['org','cov']).agg({'base':np.sum})

	bed = pd.read_csv(inputBED, sep='\t', names=['id', 'start', 'end', 'depthCoverage'])
	bed.id = bed.id.apply(lambda x: x.split("|")[1] if "|" in x else x)
	bed.start = bed.start.astype(np.int64)
	bed.end = bed.end.astype(np.int64)
	bed.depthCoverage = bed.depthCoverage.astype(np.int64)
	meta = pd.DataFrame()

	for id, ds in bed.groupby('id'):
		ds['baseCovered'] = ds['end']-ds['start']
		try:
			ds['nucleotidesCovered']=(ds['baseCovered']*ds['depthCoverage'])#/float(genomeMap[genomeMap.id==id].len)
		except:
			print inputBED, id
		meta = meta.append({'id':id, 'baseCovered' : np.sum(ds.baseCovered), 'nucleotidesCovered' : np.sum(ds.nucleotidesCovered), 'totalReads' : np.sum(ds.depthCoverage) },ignore_index=True)	

	try:
		meta = pd.merge(meta,genomeMap)
	except:
		print inputBED
		raise BaseException() 
		
	meta['depthCoverage'] = meta["totalReads"] / meta["len"]
	meta=meta.groupby('org').agg({'totalReads':np.sum, 'baseCovered':np.sum, 'depthCoverage': np.mean}).reset_index() #dividi
	meta.baseCovered=meta.baseCovered.astype(np.int64)
	meta.totalReads=meta.totalReads.astype(np.int64)
	glength = genomeMap.groupby('org', as_index=False).agg({'len':np.sum})

	out=pd.merge(glength, meta)
	out['depthCoverage']=out.totalReads/out.len
	out['breadthCoverage']=out.baseCovered/out.len
	out.to_csv(inputBED.replace(".bed","_mapped.bed"), sep='\t', index=False, float_format='%1.4f')
	relativeAbundances.calculateRelativeAbb(inputBED.replace(".bed","_mapped.bed"), contigGenome)
	
if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument("inputBED", help="BED to be mapped")
	parser.add_argument("contigGenomeCSV", help="Map contig-genome")
	args = parser.parse_args()
	if os.stat(args.inputBED).st_size > 0:
		bedmap(args.inputBED, args.contigGenomeCSV)