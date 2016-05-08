#!/usr/bin/env python 
import Bio, os, argparse, glob, sys, pickle
import pandas as pd
import multiprocessing as mp
import warnings
with warnings.catch_warnings():
	warnings.simplefilter("ignore")
	from Bio import SearchIO, SeqIO, AlignIO, SeqRecord, Seq
	from StringIO import StringIO
	from Bio.SeqRecord import SeqRecord
	from Bio.Seq import Seq
	from Bio.Blast.Applications import NcbiblastnCommandline
	from Bio.Align.Applications import MuscleCommandline
	from Bio.Phylo.Applications import RaxmlCommandline
from collections import defaultdict

LEN_THRESHOLD = 500
IDEN_THRESHOLD = 0.7
GENOME_CONSERVATION = 1
EVALUE = 1e-10
NPROC=20
BOOTSTRAP=50

muscle_exe = "/scratch/sharedCM/users/beghini/muscle3.8.1551"
raxml_exe = "/scratch/sharedCM/users/beghini/raxmlHPC-PTHREADS-SSE3"
trimal_exe = "/home/francesco.beghini/bin/trimal"

genomes = []
ref_gen = ""

def j(x):
	if(x.hit_range[1]-x.hit_range[0]<0):
		x.hit = x.hit[::-1]

def initt(terminating_):
	# This places terminating in the global namespace of the worker subprocesses.
	# This allows the worker function to access `terminating` even though it is
	# not passed as an argument to the function.
	global terminating
	terminating = terminating_

def common_seq(gene):
	if not terminating.is_set():
		query = "extractedgenes/%s.fasta" % (gene.id)
		SeqIO.write(gene, query, "fasta")
	
		records = []
		#report = defaultdict(lambda: defaultdict(dict))
		report = {}
		_tmp={}
		#print "BLASTing %s" % (gene.id)
		for genome in genomes:
			genome = os.path.split(genome)[1].split('.')[0]
			# print "\t%s..." % (genome)
			if os.path.exists("blastout/%s_%s.xml" % (genome, gene.id)):
				blast_result = SearchIO.read("blastout/%s_%s.xml" % (genome, gene.id), 'blast-xml')
 				if blast_result.param_evalue_threshold != EVALUE:
 					blastn_cline = NcbiblastnCommandline(query=query, db="blastdb/%s" % genome, evalue=EVALUE, outfmt= 5, word_size = 9, perc_identity=IDEN_THRESHOLD*100, out= "blastout/%s_%s.xml" % (genome, gene.id))
 					stdout, stderr = blastn_cline()
 			else:
 				blastn_cline = NcbiblastnCommandline(query=query, db="blastdb/%s" % genome, evalue=EVALUE, outfmt= 5, word_size = 9, perc_identity=IDEN_THRESHOLD*100, out= "blastout/%s_%s.xml" % (genome, gene.id))
 				stdout, stderr = blastn_cline()
			for blast_result in SearchIO.parse("blastout/%s_%s.xml" % (genome, gene.id), 'blast-xml'):
				if len(blast_result.hits) == 0:
					os.remove("blastout/%s_%s.xml" % (genome, gene.id))
					continue
				filtered_hits = blast_result.hsp_filter(lambda hsp : hsp.aln_span > LEN_THRESHOLD)
				if(len(filtered_hits)>0):
					filtered_hits = filtered_hits[0][0]
					filtered_hits.hit.id = os.path.split(genome)[1].split('.')[0].split('.')[0]
					filtered_hits.hit.name = filtered_hits.hit.description = ""
					#map(j,filtered_hits)
				records.extend(map(lambda HSPFragment : HSPFragment.hit, filtered_hits.fragments))
			_tmp[genome]=0
		for record in records:
			_tmp[record.id] = 1
		report[gene.id]=_tmp
		#print "%s: found in %i/%i " % (gene.id, len(records), len(genomes))
		if len(records) >= len(genomes) * GENOME_CONSERVATION:
			# print "%s is a core gene" % (gene.id)
			for record in records:
				record.seq = Bio.Seq.Seq(str(record.seq).replace('-','N'), Bio.Alphabet.DNAAlphabet)
			SeqIO.write(records, 'coregenes/%s.fasta' % (gene.id.replace(":","_")),'fasta')

		else:
			None
			# print "%s is not a core gene" % (gene.id)
	else:
		terminating.set()
	return report

def muscle_aln():
	mergedaln = {}
	# Reuse alignment
	if len(glob.glob("coregenes/*.fasta")) > 0:
		print "Running muscle on extracted sequences..."
		coregenes = glob.glob("coregenes/*.fasta")
		for gene in coregenes:
			print '\r{0:.1%} completed'.format(float(coregenes.index(gene))/len(coregenes)),
			muscle_cline = MuscleCommandline(muscle_exe, input=gene, maxiters=1)
			stdout, stderr  = muscle_cline()
			alignment = AlignIO.read(StringIO(stdout), 'fasta')
			AlignIO.write(alignment, open(gene+".aln","w"),"fasta")
			# ref_gen_len = len([x for x in list(alignment) if x.id == ref_gen][0].seq)
			for seq in alignment:
				if seq.id not in mergedaln:
					mergedaln[seq.id] = seq
				else:
					mergedaln[seq.id] += seq
			missing_genome = list(set([os.path.split(genome)[1].split('.')[0] for genome in genomes]) - set([seq.id for seq in alignment]))
			for genome in missing_genome:
				if genome not in mergedaln:
					mergedaln[genome] = SeqRecord(Seq("-" * alignment.get_alignment_length(),Bio.Alphabet.SingleLetterAlphabet()),id=genome, name=genome,description=genome)
				else:
					mergedaln[genome] += "-" * alignment.get_alignment_length()
		pickle.dump(mergedaln, open("mergedaln","wb"))			
		SeqIO.write(mergedaln.values(), "muscleout.aln", "fasta")
		trimal_cline = "%s -in muscleout.aln -out trimmed_muscleout.aln -gappyout" % (trimal_exe)
		os.system(trimal_cline)
	else:
		print "No sequences to align"
		exit(0)

def generate_phylo():
	[os.remove(x) for x in glob.glob("RAxML_*")]
	print "Generating phylogenetic tree using the multiple algnment output using RAxML..."
	try:
		#best_likelihood
		print "1/3"
		raxml_cline = RaxmlCommandline(sequences="trimmed_muscleout.aln", model="GTRGAMMA", threads=NPROC, cmd=raxml_exe, name="T1", parsimony_seed=12345, num_replicates=BOOTSTRAP)
		stdout, stderr  = raxml_cline()
		#bootstrap search
		print "2/3"
		raxml_cline = RaxmlCommandline(sequences="trimmed_muscleout.aln", model="GTRGAMMA", threads=NPROC, cmd=raxml_exe, name="T2", parsimony_seed=12345, num_replicates=BOOTSTRAP, bootstrap_seed=12345)
		stdout, stderr  = raxml_cline()
		#draw bipartition
		print "3/3"
		raxml_cline = RaxmlCommandline(sequences="trimmed_muscleout.aln", model="GTRCAT", threads=NPROC, cmd=raxml_exe, name="T3.nwk", algorithm="b", starting_tree="RAxML_bestTree.T1", bipartition_filename="RAxML_bootstrap.T2")
		stdout, stderr  = raxml_cline()
	except Exception,e: 
		print str(e)
		exit(0)

if __name__ == '__main__':
	warnings.filterwarnings("ignore")
	parser = argparse.ArgumentParser()

	parser.add_argument("-input_gff")
	parser.add_argument("-input_fna")
	parser.add_argument("-input_ffn")
	parser.add_argument("-genomes", required=True)
	parser.add_argument("-processor", default=10)
	parser.add_argument("-seq_identity", default=0.7)
	parser.add_argument("-len_threshold", default=500)
	parser.add_argument("-genome_conservation", default=0.9)
	parser.add_argument("-evalue", default=1e-50)
	parser.add_argument("-bootstrap", default=50)

	if len(sys.argv)==1:
		parser.print_help()
		sys.exit(1)
	args = parser.parse_args()

	input_gff = args.input_gff
	input_fna = args.input_fna
	input_ffn = args.input_ffn

	NPROC = int(args.processor)
	IDEN_THRESHOLD = float(args.seq_identity)
	LEN_THRESHOLD = int(args.len_threshold)
	GENOME_CONSERVATION = float(args.genome_conservation)
	EVALUE = float(args.evalue)
	BOOTSTRAP=int(args.bootstrap)

	genomes = list(set([y for x in map(glob.glob,args.genomes.split(',')) for y in x]))
	#print genomes
	for genome in genomes:
		if not os.path.exists("blastdb"):
			os.mkdir("blastdb")
		if not os.path.isfile("blastdb/"+os.path.split(genome)[1].split('.')[0]+".nhr"):
			print genome
			os.system("makeblastdb -in %s -out blastdb/%s -dbtype nucl" % (genome, os.path.split(genome)[1].split('.')[0]))
	if(input_ffn is None):
		ref_gen = os.path.split(input_fna)[1].split('.')[0]
		os.system('grep -Pv "\tmRNA\t|\tCDS\t|\texon\t|\tregion\t" %s > %s' % (input_gff, os.path.splitext(input_gff)[0]+"_genes.gff" ))	
		os.system("bedtools getfasta -fi %s -bed %s -fo %s -split"  % (input_fna, os.path.splitext(input_gff)[0]+"_genes.gff", os.path.splitext(input_fna)[0]+".genes"))
		genes = SeqIO.parse(os.path.splitext(input_fna)[0]+".genes", "fasta")
		genes = [gene for gene in genes if len(gene.seq)>500]
	else:
		genes = SeqIO.parse(input_ffn, "fasta")
		genes = [gene for gene in genes if len(gene.seq)>500]
	
	if not os.path.exists("extractedgenes"):
		os.mkdir("extractedgenes")
	if not os.path.exists("coregenes"):
		os.mkdir("coregenes")
	if not os.path.exists("blastout"):
		os.mkdir("blastout")
			
	[os.remove(x) for x in glob.glob("coregenes/*")]

	terminating = mp.Event()
	pool = mp.Pool(initializer = initt, initargs=(terminating, ), processes = NPROC)
	processes = {}
	try:
		for k,v in enumerate(pool.imap_unordered(common_seq, genes),1):
			c = processes.update(v)
			print '\r{0:.1%} completed. {1} core genes found.'.format(float(k)/len(genes), len(glob.glob("coregenes/*"))),
		print ''
		pool.close()
	except Exception, e:
		print "Failed", str(e)
		pool.terminate()
		pool.join()	
	pool.join()

	if(len(processes)>0):
		table = pd.DataFrame.from_dict(processes).transpose()
		pd.DataFrame.to_csv(table, 'gene_summary.csv')

	muscle_aln()
	if os.path.exists("muscleout.aln"):
		generate_phylo()
		print "Done. The final phylogenetic tree is RAxML_bipartitionsBranchLabels.T3.nwk"
	else:
		print "Muscle alignment failed"
		exit(0)
