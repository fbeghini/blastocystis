
exclude=("Sample_752-r1.fasta.gz" "Sample_781-r1.fasta.gz" "Sample_785-r1.fasta.gz" "Sample_792-r1.fasta.gz" "Sample_804-r2.fasta.gz" "Sample_804-r1.fasta.gz")
samples=`find "/banche_dati/sharedCM/Dropbox (CIBIO)/mycobacteria/data/genomes/" -type f | grep ".fasta.gz" | sort`

echo $samples

for ef in ${exclude[*]} ;
do
  echo ${ef}
  echo $sample
  samples=`echo ${samples} | grep -v ${ef}`
  #$sagrep -v $ef "$samples"
done
echo $samples
