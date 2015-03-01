exclude=("Sample_752-r1.fasta.gz" "Sample_781-r1.fasta.gz" "Sample_785-r1.fasta.gz" "Sample_792-r1.fasta.gz" "Sample_804-r2.fasta.gz" "Sample_804-r1.fasta.gz")
samples=$(find "/banche_dati/sharedCM/Dropbox (CIBIO)/mycobacteria/data/genomes/" -type f | grep ".fasta.gz" | sort)

#for ef in ${exclude[*]}; do
#  echo $ef
  #samples=`echo ${samples} | grep -v ${ef}`
#  $sagrep -v $ef "$samples"
#done
#echo $samples

#touch /home/francesco.beghini/mycobacteria/mycobacteria_samples.fa
#cd "/banche_dati/sharedCM/Dropbox (CIBIO)/mycobacteria/data/genomes"
#for sample in $(cat ~/samples.txt)
#do
#  echo $sample
  #zcat $sample >> /home/francesco.beghini/mycobacteria/mycobacteria_samples.fa
#donea

#du /home/francesco.beghini/mycobacteria/mycobacteria_samples.fa

cd "/home/francesco.beghini"
# Build bowtie2 index and inspect

#bowtie2-build mycobacteria/mycobacteria_samples.fa mycobacteria/mycobacteria
#samtools faidx mycobacteria/mycobacteria_samples.fa

metagenomes=("G18282" "G18107" "G18131" "G18221" "G18246" "G18180" "G18299")

for mg in ${metagenomes[*]}; do
  #mkdir $mg
  #tar xjf /banche_dati/sharedCM/metagenomes/hmpii/$mg.tar.bz2 -O | \
  #pyphlan/fastx_len_filter.py --min_len 80 | \
  #bowtie2 -a --no-unal --very-sensitive -x mycobacteria/mycobacteria -U - -p 4 | \
  #samtools view -Sb - > $mg/$mg.bam
  #samtools sort $mg/$mg.bam $mg/$mg.srtd.bam -m 12000000000
  rm $mg/$mg.bam
  mv $mg/$mg.srtd.bam $mg/$mg.bam
  #samtools index $mg/$mg.bam 
  samtools mpileup -uf mycobacteria/mycobacteria_samples.fa $mg/$mg.bam | \
  bcftools view -bvcg - > $mg/$mg.bcf
  genomeCoverageBed  -bg -max 1 -ibam $mg/$mg.bam -g mycobacteria/mycobacteria_samples.fa.fai > $mg/$mg_coverage.tsv
done
