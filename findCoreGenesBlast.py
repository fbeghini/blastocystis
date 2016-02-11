import Bio, os, argparse, glob
from Bio import SearchIO, SeqIO, AlignIO
from StringIO import StringIO
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Align.Applications import MuscleCommandline
from Bio.Phylo.Applications import RaxmlCommandline
from collections import defaultdict
import multiprocessing as mp
import pickle

LEN_THRESHOLD = 500
IDEN_THRESHOLD = 0.7
muscle_exe = "/scratch/sharedCM/users/beghini/muscle3.8.1551"
raxml_exe = "/scratch/sharedCM/users/beghini/raxmlHPC-PTHREADS-SSE3"

report = defaultdict(lambda: defaultdict(dict))

def common_seq(gene, genomes):
	query = "extractedgenes/%s.fasta" % (gene.id)
	if not os.path.exists("extractedgenes"):
		os.mkdir("extractedgenes")
	if not os.path.exists("coregenes"):
		os.mkdir("coregenes")
	if not os.path.exists("blastout"):
		os.mkdir("blastout")
	SeqIO.write(gene, query, "fasta")
	records = []
	print "BLASTing %s on " % (gene.id)
	for genome in genomes:
		print "\t%s..." % (genome)
		blastn_cline = NcbiblastnCommandline(query=query, db="blastdb/%s" % genome.split('.')[0], evalue=0.001, outfmt= 5, word_size = 9, out= "blastout/%s.xml" % gene.id)# % genome.split('.')[0])
		stdout, stderr = blastn_cline()
		if len(stderr) == 0:										# handle this
			for blast_result in SearchIO.parse("blastout/%s.xml" % gene.id, 'blast-xml'):
				filtered_hits = blast_result.hsp_filter(lambda hsp : hsp.aln_span > LEN_THRESHOLD)
				filtered_hits = filtered_hits.hsp_filter(lambda hsp : hsp.ident_num/float(hsp.aln_span) > IDEN_THRESHOLD)
				if(len(filtered_hits)>0):
					filtered_hits = filtered_hits[0][0]
					filtered_hits.hit.id = genome.split('.')[0]
					filtered_hits.hit.name = filtered_hits.hit.description = ""
				records.extend(map(lambda HSPFragment : HSPFragment.hit, filtered_hits.fragments))
	for record in records:
		report[gene.id][record.id] = 1

	if len(records) >= len(genomes):
		print "%s is a core gene" % (gene.id)
		for record in records:
			record.seq = Bio.Seq.Seq(str(record.seq).replace('-','N'), Bio.Alphabet.DNAAlphabet)
		SeqIO.write(records, 'coregenes/%s.fasta' % (gene.id.replace(":","_")),'fasta')

	else:
		print "%s is not a core gene" % (gene.id)

def muscle_aln():
	# if not os.path.exists("musclealn"):
		# os.mkdir("musclealn")
	mergedaln = {}
	print "Running muscle on extracted sequences..."
	for gene in glob.glob("coregenes/*.fasta"):
		muscle_cline = MuscleCommandline(muscle_exe, input=gene) 
		stdout, stderr  = muscle_cline()
		alignment = AlignIO.read(StringIO(stdout), 'fasta')
		for seq in alignment:
			if seq.id not in mergedaln:
				mergedaln[seq.id] = seq
			else:
				mergedaln[seq.id] += seq
	SeqIO.write(mergedaln.values(), "muscleout.fasta", "fasta")

def generate_phylo():
	#os.system("google-chrome-stable http://bit.ly/1nxxSEW")
	print "Generating phylogenetic tree using the multiple algnment output using RAxML..."
	try:
		#best_likelihood
		raxml_cline = RaxmlCommandline(sequences="muscleout.fasta", model="GTRGAMMA", threads=20, cmd=raxml_exe, name="T1", parsimony_seed=12345, num_replicates=5)
		stdout, stderr  = raxml_cline()
		#bootstrap search
		raxml_cline = RaxmlCommandline(sequences="muscleout.fasta", model="GTRGAMMA", threads=20, cmd=raxml_exe, name="T2", parsimony_seed=12345, num_replicates=5, bootstrap_seed=12345)
		stdout, stderr  = raxml_cline()
		#draw bipartition
		raxml_cline = RaxmlCommandline(sequences="muscleout.fasta", model="GTRCAT", threads=20, cmd=raxml_exe, name="T3.nwk", algorithm="b", starting_tree="RAxML_bestTree.T1", bipartition_filename="RAxML_bootstrap.T2")
		stdout, stderr  = raxml_cline()
	except :
		None

if __name__ == '__main__':
	parser = argparse.ArgumentParser()

	parser.add_argument("-input_gff")
	parser.add_argument("-input_fna")
	parser.add_argument("-genomes")
	args = parser.parse_args()

	input_gff = args.input_gff #'/CIBIO/sharedCM/projects/blastocystis/dataset/Blastocystis/ST4/GCF_000743755.1_ASM74375v1_genomic.gff'
	input_fna = args.input_fna #'/CIBIO/sharedCM/projects/blastocystis/dataset/Blastocystis/ST4/GCF_000743755.1_ASM74375v1_genomic.fna'
	genomes = args.genomes.split(',') #["GCF_000151665.1_ASM15166v1_genomic.fasta", "GCF_000743755.1_ASM74375v1_genomic.fasta", "ST2_JZRJ01.1.fasta", "ST3_JZRK01.1.fasta", "ST6_JZRM01.1.fasta", "ST8_JZRN01.1.fasta", "ST9_JZRO01.1.fasta"]
	
	for genome in genomes:
		if not os.path.exists("blastdb"):
			os.mkdir("blastdb")
		if not os.path.isfile("blastdb/"+genome.split('.')[0]+".nhr"):
			os.system("makeblastdb -in %s -out blastdb/%s -dbtype nucl" % (genome, genome.split('.')[0]))

	os.system('grep -Pv "\tmRNA\t|\tCDS\t|\texon\t|\tregion\t" %s > %s' % (input_gff, os.path.splitext(input_gff)[0]+"_genes.gff" ))	
	os.system("bedtools getfasta -fi %s -bed %s -fo %s -split"  % (input_fna, os.path.splitext(input_gff)[0]+"_genes.gff", os.path.splitext(input_fna)[0]+"_genes.fna"))

	genes = SeqIO.parse(os.path.splitext(input_fna)[0]+"_genes.fna", "fasta")

	pool = mp.Pool(processes = 50)
	processes = [pool.apply(common_seq, args=(gene.upper(), genomes)) for gene in genes]
	pool.close()
	pool.terminate()
	pickle.dump(report,"report.pkl")
	muscle_aln()
	generate_phylo()
	print "Done. The final phylogenetic tree is RAxML_bipartitionsBranchLabels.T3.nwk"