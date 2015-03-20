#!/usr/bin/env python

import os, argparse, sys, subprocess, tarfile
from Bio import SeqIO
from tempfile import NamedTemporaryFile

parser = argparse.ArgumentParser()
parser.add_argument("--ref", help="FASTA containing the reference genoma")
parser.add_argument("--input-type", help="FR for the Illumina paired ends, U for unpaired reads, SRA", choices=["FR","U", "SRA"])
parser.add_argument("--metagenomes", help="Comma-separated list of metagenomes to analyze")
parser.add_argument("--basename-index", help="Basename for the index")
parser.add_argument("--output-folder", help="Output folder")
parser.add_argument("--force", help="Force the rebuild of the index", action='store_true')
args = parser.parse_args()

tempdir = "/scratch/sharedCM/users/beghini"

if args.ref:
	#Build bowtie2 index and inspect
	if args.force:
		print "Building the index..."
		subprocess.Popen("bowtie2-build %s %s" % (args.ref, args.basename_index), shell=True)
		subprocess.Popen("samtools faidx %s" % (args.ref))

	for mg in args.metagenomes.split(','):
		outname = args.output_folder + "/" + mg.split('/')[-1].split('.')[0]
		#print outname	
		with NamedTemporaryFile(delete=True) as _bam:
			#Extract mates for paired inputs with FR
			if(args.input_type == 'FR'):
				print "Listing the content of %s ..." %mg
				mates = tarfile.open(mg, 'r:bz2').getnames()
				mate1 = [m for m in mates if "_1" in m]
				mate2 = [m for m in mates if "_2" in m]

				forward = NamedTemporaryFile(delete=True, dir=tempdir)
				reverse = NamedTemporaryFile(delete=True, dir=tempdir)

				print "\tBuilding forward mate"
				tf = subprocess.Popen("tar -xjf %s %s -O"%(mg, " ".join(mate1)), shell=True, stdout=forward)
				tf.wait()
				print "\tBuiling reverse mate"
				tr = subprocess.Popen("tar -xjf %s %s -O"%(mg, " ".join(mate2)), shell=True, stdout=reverse)
				tr.wait()
				print "Done."
				com = "bowtie2 --no-unal -a --very-sensitive -p 6 -x %s -1 %s -2 %s | samtools view -Sb -" % (args.basename_index, forward.name, reverse.name)
			# bzip the unpaired input and use it as stdin for bowtie2
			elif(args.input_type == 'U'):
				com = "tar -xjv %s -O | bowtie2 --no-unal -a --very-sensitive -p 6 -x %s -U - | samtools view -Sb -" % (mg,args.basename_index)
			elif(args.input_type == 'SRA'):
				com = "fastq-dump -Z %s --split-spot | bowtie2 --no-unal -a --very-sensitive -p 6 -x %s -U - | samtools view -Sb -" % (mg,args.basename_index)
			#print com
			print "Aligning and sorting %s ..." %mg
			bam = subprocess.Popen(com, shell=True, stdout=_bam)
			bam.wait()
			sort =subprocess.Popen("samtools sort -	 %s -@ 4 -m 6G" % (outname), shell=True, stdin=_bam)
			sort.wait()
			if args.input_type=="FR":
				forward.close()
				reverse.close()
		#print outname

		index = subprocess.Popen("samtools index %s.bam" % (outname), shell=True)
		index.wait()
		
		# print "samtools sort - %s.bam" % (outname)
		# print "samtools index %s.bam" % (outname)
		# print "samtools mpileup -uf %s %s.bam | bcftools view -bvcg -" % (args.ref, outname)
		
		mpileup = subprocess.Popen("samtools mpileup -uf %s %s.bam | bcftools view -bvcg -" % (args.ref, outname), shell=True, stdout=subprocess.PIPE)
	
		with open("%s.bcf" % (outname),"w") as bcfout:
			bcfout.writelines(mpileup.stdout)

		bed = subprocess.Popen("bedtools genomecov -ibam %s.bam -g %s.fai" % (outname, args.ref), shell=True, stdout=subprocess.PIPE)
		bed.wait()

		g = subprocess.Popen('grep -P "\t0\t"', shell=True, stdin=bed.stdout)
		
		with open("%s.tsv" % (outname),"w") as tsvout:
			tsvout.writelines(g.stdout)
		
		bed = subprocess.Popen("bedtools genomecov -bg -ibam %s.bam -g %s.fai" % (outname, args.ref), shell=True, stdout=subprocess.PIPE)
		
		with open("%s.bed" % (outname),"w") as bedout:
			bedout.writelines(bed.stdout)
		# print "bedtools genomecov -bg -max 1 -ibam %s.bam -g %s.fai" % (outname, args.ref)