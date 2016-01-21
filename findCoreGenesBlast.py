import Bio
from Bio import SearchIO, SeqIO
from Bio.Blast.Applications import NcbiblastnCommandline
import os, argparse
import multiprocessing as mp

def common_seq(gene, genomes):
	query = "extractedgenes/%s.fasta" % (gene.id)
	if not os.path.exists("extractedgenes"):
		os.mkdir("extractedgenes")
	if not os.path.exists("coregenes"):
		os.mkdir("coregenes")
	SeqIO.write(gene, query, "fasta")
	records = []
	print "BLASTing %s on " % (gene.id)
	for genome in genomes:
		print "\t%s..." % (genome)
		blastn_cline = NcbiblastnCommandline(query=query, db="blastdb/%s" % genome.split('.')[0], evalue=0.001, outfmt= 5, word_size = 9, out= "%s.xml" % gene.id)# % genome.split('.')[0])
		stdout, stderr = blastn_cline()
		for blast_result in SearchIO.parse("%s.xml" % gene.id, 'blast-xml'):
			filtered_hits = blast_result.hsp_filter(lambda hsp : hsp.aln_span > LEN_THRESHOLD)
			filtered_hits = filtered_hits.hsp_filter(lambda hsp : hsp.ident_num/float(hsp.aln_span) > IDEN_THRESHOLD)
			records.extend(map(lambda HSPFragment : HSPFragment.hit, filtered_hits.fragments))
	if len(records) >= len(genomes):
		print "%s is a core gene" % (gene.id)
		for record in records:
			record.seq = Bio.Seq.Seq(str(record.seq).replace('-','N'), Bio.Alphabet.DNAAlphabet)
		SeqIO.write(records, 'coregenes/%s.fasta' % (gene.id),'fasta')
	else:
		print "%s is not a core gene" % (gene.id)
	os.remove("%s.xml" % gene.id)

if __name__ == '__main__':
	LEN_THRESHOLD = 500
	IDEN_THRESHOLD = 0.7

	input_gff = '/CIBIO/sharedCM/projects/blastocystis/dataset/Blastocystis/ST4/GCF_000743755.1_ASM74375v1_genomic.gff'
	input_fna = '/CIBIO/sharedCM/projects/blastocystis/dataset/Blastocystis/ST4/GCF_000743755.1_ASM74375v1_genomic.fna'
	genomes = ["GCF_000151665.1_ASM15166v1_genomic.fasta", "GCF_000743755.1_ASM74375v1_genomic.fasta", "ST2_JZRJ01.1.fasta", "ST3_JZRK01.1.fasta", "ST6_JZRM01.1.fasta", "ST8_JZRN01.1.fasta", "ST9_JZRO01.1.fasta"]
	for genome in genomes:
		if not os.path.exists("blastdb"):
			os.mkdir("blastdb")
		if not os.path.isfile("blastdb/"+genome.split('.')[0]+".nhr"):
			os.system("makeblastdb -in %s -out blastdb/%s -dbtype nucl" % (genome, genome.split('.')[0]))

	os.system('grep -Pv "\tmRNA\t|\tCDS\t|\texon\t|\tregion\t" %s > %s' % (input_gff, os.path.splitext(input_gff)[0]+"_genes.gff" ))
	os.system("bedtools getfasta -fi %s -bed %s -fo %s -split"  % (input_fna, os.path.splitext(input_gff)[0]+"_genes.gff", os.path.splitext(input_fna)[0]+"_genes.fna"))

	genes = SeqIO.parse(os.path.splitext(input_fna)[0]+"_genes.fna", "fasta")

	pool = mp.Pool(processes=10)
	processes = [pool.apply(common_seq, args=(gene.upper(), genomes)) for gene in genes]