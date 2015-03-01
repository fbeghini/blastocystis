import os, argparse, sys, subprocess, tarfile
from Bio import SeqIO
from tempfile import NamedTemporaryFile

parser = argparse.ArgumentParser()
parser.add_argument("--ref", help="FASTA containing the reference genoma")
parser.add_argument("--input-type", help="FR for the Illumina paired ends, U for unpaired reads", choices=["FR","U"])
parser.add_argument("--metagenomes", help="Comma-separated list of metagenomes to analyze")
parser.add_argument("--basename-index", help="Basename for the index")
parser.add_argument("--output-folder", help="Output folder")
parser.add_argument("--force", help="Force the rebuild of the index", action='store_true')
args = parser.parse_args()

#def buildDataset(sequences):
	# try:
	# 	for s in sequences:
	# 		subprocess.Popen("args")
if args.ref:
	#Build fasta dataset

	#Build bowtie2 index and inspect
	if args.force:
		print "Building the index..."
		subprocess.Popen("bowtie2-build %s %s" % (args.ref, args.basename_index), shell=True)
		subprocess.Popen("samtools faidx %s" % (args.ref))

	for mg in args.metagenomes.split(','):
		outname = args.output_folder + "/" + mg.split('/')[-1].split('.')[0]
		# print outname
		#Extract mates for paired inputs with FR
		with NamedTemporaryFile(delete=True) as _bam:
			if(args.input_type == 'FR'):
				print "Listing the content of %s ..." %mg
				mates = tarfile.open(mg, 'r:bz2').getnames()
				mate1 = [m for m in mates if "_1" in m]
				mate2 = [m for m in mates if "_2" in m]

				print "Alignment..."
				bam = subprocess.Popen(
					"bowtie2 --no-unal -a --very-sensitive -p 4 -x %s "
					"-1 <(tar -xjf %s {%s} -O) "
					"-2 <(tar -xjf %s {%s} -O) "
					"| samtools view -Sb -" % (args.basename_index,mg, ','.join(mate1), mg, ','.join(mate2))
					, shell=True, executable='/bin/bash', stdout=_bam)
				# print "bowtie2 --no-unal -a --very-sensitive -p 4 -x %s/%s -1 <(tar -xjf %s {%s} -O) -2 <(tar -xjf %s {%s} -O) | samtools view -Sb -" % (args.basename_index,args.basename_index,mg, ",mate1", mg, ",mate2")
			# bzip the unpaired input and use it as stdin for bowtie2
			elif(args.input_type == 'U'):
				print "Alignment..."
				bam = subprocess.Popen(
					"tar -xjv %s -O | "
					"bowtie2 --no-unal -a --very-sensitive -p 4 -x %s -U - | "
					"samtools view -Sb -" % (mg,args.basename_index)
					, shell=True, executable='/bin/bash', stdout=_bam)
			bam.wait()
			subprocess.Popen("samtools sort - %s -m 12000000000" % (outname), shell=True, stdin=_bam)

		subprocess.Popen("samtools index %s.bam" % (outname), shell=True)
		
		# print "samtools sort - %s.bam" % (outname)
		# print "samtools index %s.bam" % (outname)
		# print "samtools mpileup -uf %s %s.bam | bcftools view -bvcg -" % (args.ref, outname)
		
		mpileup = subprocess.Popen("samtools mpileup -uf %s %s.bam | bcftools view -bvcg -" % (args.ref, outname), shell=True, stdout=subprocess.PIPE)
		with open("%s.bcf" % (outname),"w") as bcfout:
			bcfout.writelines(mpileup.stdout)

		bed = subprocess.Popen("bedtools genomecov -ibam %s.bam -g %s.fai" % (outname, args.ref), shell=True, stdout=subprocess.PIPE)
		with open("%s_coverage.tsv" % (outname),"w") as tsvout:
			tsvout.writelines(bed.stdout)

		# print "bedtools genomecov -bg -max 1 -ibam %s.bam -g %s.fai" % (outname, args.ref)