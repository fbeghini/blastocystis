#! /usr/bin/env python

import os, argparse, sys, subprocess, tarfile, fastq_len_filter, StringIO, map_bed
from Bio import SeqIO
from tempfile import NamedTemporaryFile

bowtie2_path = "/CIBIO/sharedCM/projects/mycobacteria/bowtie2-2.2.5/"

def clean_files(input_type, tempdir, mgname):
	if input_type=="SRA":
		os.remove(tempdir+"/"+mgname+"_1.fastq")
		os.remove(tempdir+"/"+mgname+"_2.fastq")
		os.remove(tempdir+"/"+mgname+"_1.filtered.fastq")
		os.remove(tempdir+"/"+mgname+"_2.filtered.fastq")
	if input_type=="FR":
		forward.close()
		reverse.close()

parser = argparse.ArgumentParser()

parser.add_argument("--ref", help="FASTA containing the reference genoma", required=True)
parser.add_argument("--input-type", help="FR for paired ends, U for unpaired reads, SRA", choices=["FR","U", "SRA"], required=True)
parser.add_argument("--metagenomes", help="Comma-separated list of metagenomes to analyze", required=True)
parser.add_argument("--basename-index", help="Basename for the index", required=True)
parser.add_argument("--len-filter", type = int, default=90)
parser.add_argument("--genomeMap", required=True)
parser.add_argument("--output-folder", help="Output folder", required=True)

args = parser.parse_args()
tempdir = "/scratch/sharedCM/users/beghini"
#if args.ref and args.input_type and args.metagenomes and args.basename_index :

mg = args.metagenomes
outname = os.path.abspath(args.output_folder) + "/" + mg.split('/')[-1].split('.')[0]
#Extract mates for paired inputs with FR

if(args.input_type == 'FR'):
	print "Listing the content of %s ..." %mg
	try:
		mates = tarfile.open(mg, 'r:bz2').getnames()
	except IOError, e:
		raise e
	mate1 = [m for m in mates if "_1" in m]
	mate2 = [m for m in mates if "_2" in m]

	forward	 = NamedTemporaryFile(delete=True, dir=tempdir)
	reverse = NamedTemporaryFile(delete=True, dir=tempdir)

	print "\tBuilding forward mate..."
	tf = subprocess.Popen("tar -xjf %s %s -O"%(mg, " ".join(mate1)), shell=True, stdout=subprocess.PIPE)
	fastq_len_filter.filter(args.len_filter, StringIO.StringIO(tf.communicate()[0]),forward)

	print "\tBuiling reverse mate..."
	tr = subprocess.Popen("tar -xjf %s %s -O"%(mg, " ".join(mate2)), shell=True, stdout=subprocess.PIPE)
	fastq_len_filter.filter(args.len_filter, StringIO.StringIO(tr.communicate()[0]),reverse)

	print "Done!"
	com = "bowtie2 --no-unal -a --very-sensitive -p 2 -x %s -1 %s -2 %s | samtools view -Sb -" % (args.basename_index, forward.name, reverse.name)

elif(args.input_type == 'U'):
	com = "bowtie2 --no-unal -a --very-sensitive -p 2 -x %s -U - | samtools view -Sb -" % (args.basename_index)

elif(args.input_type == 'SRA'):
	print "fastq dump of %s..." % mg
	dump = subprocess.Popen("fastq-dump %s --split-3 -O %s" % (mg, tempdir), shell=True)
	dump.wait()
	mgname = mg.split(".")[0]

	fastq_len_filter.filter(args.len_filter,open(tempdir+"/"+mgname+"_1.fastq",'r'),open(tempdir+"/"+mgname+"_1.filtered.fastq",'w'))
	fastq_len_filter.filter(args.len_filter,open(tempdir+"/"+mgname+"_2.fastq",'r'),open(tempdir+"/"+mgname+"_2.filtered.fastq",'w'))
	
	com = "bowtie2 --no-unal -a --very-sensitive -p 2 -x %s -1 %s_1.filtered.fastq -2 %s_2.filtered.fastq | samtools view -Sb -" % (args.basename_index, tempdir+"/"+mgname, tempdir+"/"+mgname)

try:
	with NamedTemporaryFile(delete=True, dir=tempdir) as _bam:
		print "Aligning %s..." %mg
		bam = subprocess.Popen( bowtie2_path+com, shell=True, stdout=_bam)
		if(args.input_type=='U'):
			bam.stdin = sys.stdin
		bam.communicate()
		# bt2ret = bam.returncode
		print "Sorting %s ..." %mg
		# if bt2ret!=0: raise BaseException
		sort =subprocess.Popen("samtools sort - %s -@ 4 -m 6G" % (outname), shell=True, stdin=_bam)
		sort.communicate()
except:
	if args.input_type=="SRA":
		os.remove(tempdir+"/"+mgname+"_1.fastq")
		os.remove(tempdir+"/"+mgname+"_2.fastq")
		os.remove(tempdir+"/"+mgname+"_1.filtered.fastq")
		os.remove(tempdir+"/"+mgname+"_2.filtered.fastq")
	if args.input_type=="FR":
		forward.close()
		reverse.close()
	raise BaseException("bowtie2 failed to align %s" % mg)

if args.input_type=="SRA":
	os.remove(tempdir+"/"+mgname+"_1.fastq")
	os.remove(tempdir+"/"+mgname+"_2.fastq")
	os.remove(tempdir+"/"+mgname+"_1.filtered.fastq")
	os.remove(tempdir+"/"+mgname+"_2.filtered.fastq")
if args.input_type=="FR":
	forward.close()
	reverse.close()

try:
	index = subprocess.Popen("samtools index %s.bam" % (outname), shell=True)
	index.wait()
except:
	raise BaseException("Sorting of %s has failed" % mg)

mpileup = subprocess.Popen("samtools mpileup -uf %s %s.bam | bcftools view -bvcg -" % (args.ref, outname), shell=True, stdout=subprocess.PIPE)

with open("%s.bcf" % (outname),"w") as bcfout:
	bcfout.writelines(mpileup.communicate()[0])

bed = subprocess.Popen("bedtools genomecov -ibam %s.bam -g %s.fai" % (outname, args.ref), shell=True, stdout=subprocess.PIPE)

with open("%s.tsv" % (outname),"w") as tsvout:
	tsvout.writelines(bed.communicate()[0])

bed = subprocess.Popen("bedtools genomecov -bg -ibam %s.bam -g %s.fai" % (outname, args.ref), shell=True, stdout=subprocess.PIPE)

with open("%s.bed" % (outname),"w") as bedout:
	bedout.writelines(bed.communicate()[0])

map_bed.bedmap("%s.bed" % outname, args.genomeMap)
