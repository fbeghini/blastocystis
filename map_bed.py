#!/usr/bin/env python
import argparse, os
import numpy as np
import pandas as pd
import relativeAbundances

def bedmap(inputBED, contigGenome, contigsRemove):
	genomeMap = pd.read_csv(contigGenome, sep=';', names=['id','org','len'])
	genomeMap.id = genomeMap.id.apply(lambda x: x.split("|")[1] if "|" in x else x)
	genomeMap=genomeMap.drop_duplicates('id')
	#genomeMap.to_csv('/tmp/pangenome.csv', sep=';', header=False, index=False)
	
	# tsv = pd.read_csv(args.inputTSV, sep=';', names=['id','org','len'])
	# tsv = tsv[tsv['id']!='genome']
	# tsv.id = tsv.id.apply(lambda x: x.split("|")[1] if "|" in x else x)
	
	# contigs = pd.merge(tsv, genomeMap, how='inner', on='id').groupby(['org','cov']).agg({'base':np.sum})
        
        contigs_to_remove = pd.read_csv(contigsRemove, names=['id'])

	bed = pd.read_csv(inputBED, sep='\t', names=['id', 'start', 'end', 'depthCoverage'])
	bed.id = bed.id.apply(lambda x: x.split("|")[1] if "|" in x else x)
	bed.start = bed.start.astype(np.int64)
	bed.end = bed.end.astype(np.int64)
	bed.depthCoverage = bed.depthCoverage.astype(np.int64)
        bed = bed[~bed['id'].isin(contigs_to_remove['id'])]
	genomeMap = genomeMap[~genomeMap['id'].isin(contigs_to_remove['id'])]
        meta = pd.DataFrame()

	for id, ds in bed.groupby('id'):
		ds.loc[:,'baseCovered'] = ds['end']-ds['start']
		try:
			ds.loc[:,'nucleotidesCovered']=(ds['baseCovered']*ds['depthCoverage'])#/float(genomeMap[genomeMap.id==id].len)
		except:
			print inputBED, id
		meta = meta.append({'id':id, 'baseCovered' : np.sum(ds.baseCovered), 'nucleotidesCovered' : np.sum(ds.nucleotidesCovered), 'totalReads' : np.sum(ds.depthCoverage) },ignore_index=True)	

	try:
		meta = pd.merge(meta,genomeMap)
	except:
                meta = genomeMap.groupby('org', as_index=False).agg({'len':np.sum})
                meta['baseCovered'] = 0
                meta['totalReads'] = 0
                meta['depthCoverage'] = 0
                meta['breadthCoverage'] = 0
	        meta.to_csv(inputBED.replace(".bed","_mapped.bed"), sep='\t', index=False, float_format='%1.4f')

	meta['depthCoverage'] = meta["totalReads"] / meta["len"]
	meta=meta.groupby('org').agg({'totalReads':np.sum, 'baseCovered':np.sum, 'depthCoverage': np.mean}).reset_index() #dividi
	meta.baseCovered=meta.baseCovered.astype(np.int64)
	meta.totalReads=meta.totalReads.astype(np.int64)
	glength = genomeMap.groupby('org', as_index=False).agg({'len':np.sum})

	out=pd.merge(glength, meta, how='outer')
	out['depthCoverage']=out.totalReads/out.len
	out['breadthCoverage']=out.baseCovered/out.len
	out.to_csv(inputBED.replace(".bed","_mapped.bed"), sep='\t', index=False, float_format='%1.4f')
	#relativeAbundances.calculateRelativeAbb(inputBED.replace(".bed","_mapped.bed"), contigGenome)
	
if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument("inputBED", help="BED to be mapped")
	parser.add_argument("contigGenomeCSV", help="Map contig-genome")
	parser.add_argument("filterContigs", help="Contigs to be removed")
        args = parser.parse_args()
	if os.stat(args.inputBED).st_size > 0:
		bedmap(args.inputBED, args.contigGenomeCSV, args.filterContigs)
