#! /usr/bin/env python

import os, argparse, sys, subprocess, tarfile
from Bio import SeqIO
from tempfile import NamedTemporaryFile

parser = argparse.ArgumentParser()
parser.add_argument("--ref", help="FASTA containing the reference genoma")
parser.add_argument("--input-type", help="FR for the Illumina paired ends, U for unpaired reads, SRA, GZ", choices=["FR","U","GZ", "SRA"])
parser.add_argument("--metagenomes", help="Comma-separated list of metagenomes to analyze")
parser.add_argument("--basename-index", help="Basename for the index")
parser.add_argument("--output-folder", help="Output folder")
args = parser.parse_args()

tempdir = "/scratch/sharedCM/users/beghini"

if args.ref and args.input_type and args.metagenomes and args.basename_index :

	# if args.force:
	# 	print "Building the index..."
	# 	subprocess.Popen("bowtie2-build %s %s" % (args.ref, args.basename_index), shell=True)
	# 	subprocess.Popen("samtools faidx %s" % (args.ref))

	for mg in args.metagenomes.split(','):
		outname = os.path.abspath(args.output_folder) + "/" + mg.split('/')[-1].split('.')[0]
		#Extract mates for paired inputs with FR

		if(args.input_type == 'FR'):
			print "Listing the content of %s ..." %mg
			mates = tarfile.open(mg, 'r:bz2').getnames()
			mate1 = [m for m in mates if "_1" in m]
			mate2 = [m for m in mates if "_2" in m]

			forward = NamedTemporaryFile(delete=True, dir=tempdir)
			reverse = NamedTemporaryFile(delete=True, dir=tempdir)

			print "\tBuilding forward mate..."
			tf = subprocess.Popen("tar -xjf %s %s -O"%(mg, " ".join(mate1)), shell=True, stdout=forward)
			tf.communicate()
			print "\tBuiling reverse mate..."
			tr = subprocess.Popen("tar -xjf %s %s -O"%(mg, " ".join(mate2)), shell=True, stdout=reverse)
			tr.communicate()
			print "Done!"
			com = "bowtie2 --no-unal -a --very-sensitive -p 2 -x %s -1 %s -2 %s | samtools view -Sb -" % (args.basename_index, forward.name, reverse.name)

		elif(args.input_type == 'U'):
			com = "tar -xjf %s -O | bowtie2 --no-unal -a --very-sensitive -p 2 -x %s -U - | samtools view -Sb -" % (mg,args.basename_index)
		
		elif(args.input_type=='GZ'):
			com = "zcat %s | bowtie2 --no-unal -a --very-sensitive -p 2 -x %s -U - | samtools view -Sb -" % (mg,args.basename_index)

		elif(args.input_type == 'SRA'):
			print "fastq dump of %s..." % mg
			dump = subprocess.Popen("fastq-dump %s --split-3 -O %s" % (mg, tempdir), shell=True)
			dump.wait()
			mgname = mg.split(".")[0]
			com = "bowtie2 --no-unal -a --very-sensitive -p 2 -x %s -1 %s_1.fastq -2 %s_2.fastq | samtools view -Sb -" % (args.basename_index, tempdir+"/"+mgname, tempdir+"/"+mgname)

		with NamedTemporaryFile(delete=True, dir=tempdir) as _bam:
			bam = subprocess.Popen(com, shell=True, stdout=_bam)
			print "Aligning %s..." %mg
			bam.communicate()
			bt2ret = bam.returncode
			sort =subprocess.Popen("samtools sort - %s -@ 4 -m 6G" % (outname), shell=True, stdin=_bam)
			print "Sorting %s ..." %mg
			sort.communicate()
		
		if args.input_type=="SRA":
			os.remove(tempdir+"/"+mgname+"_1.fastq")
			os.remove(tempdir+"/"+mgname+"_2.fastq")
		
		if args.input_type=="FR":
			forward.close()
			reverse.close()

		if bt2ret==0:
			index = subprocess.Popen("samtools index %s.bam" % (outname), shell=True)
			index.wait()

			mpileup = subprocess.Popen("samtools mpileup -uf %s %s.bam | bcftools view -bvcg -" % (args.ref, outname), shell=True, stdout=subprocess.PIPE)
			
			with open("%s.bcf" % (outname),"w") as bcfout:
				bcfout.writelines(mpileup.communicate()[0])

			bed = subprocess.Popen("bedtools genomecov -ibam %s.bam -g %s.fai | sort" % (outname, args.ref), shell=True, stdout=subprocess.PIPE)
			
			with open("%s.tsv" % (outname),"w") as tsvout:
				tsvout.writelines(bed.communicate()[0])
			
			bed = subprocess.Popen("bedtools genomecov -bg -ibam %s.bam -g %s.fai | sort" % (outname, args.ref), shell=True, stdout=subprocess.PIPE)

			with open("%s.bed" % (outname),"w") as bedout:
				bedout.writelines(bed.communicate()[0])
		else:
			print "bowtie2 failed to align %s" % mg