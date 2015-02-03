#mpile=$(cat "/home/francesco.beghini/G18282/G18282_coverage.tsv")
#map="/home/francesco.beghini/contig.csv"

#while read line;
#do
	#gread=$(echo $line | awk '{print $1}')
	#genome=$(awk -F "," -v gread="$gread" '{if($1==gread) print $2}' $map)
	#echo ${mpile//$gread/$genome} >> /home/francesco.beghini/G18282/G18282_coverage.tsv.map
#done < $map


import os, csv
from collections import defaultdict

mapped=open("G18282/mapped_G18282.tsv.new","w+")

map = {}																		# Dict of {read:genome}
mappedcontigs = defaultdict(dict)
contigs = defaultdict(dict)														# Dict of dict of list of read of a genome
																				# k:v -> {contig as string : 
																				#			{coverage as int : [contig as int, total as int] }
																				#		 }
with open("contig.csv","r") as csvfile:
	reader = csv.reader(csvfile, delimiter=';')
	map = dict(reader)

with open("G18282/G18282_coverage.tsv","r") as tsv:
	reader = csv.reader(tsv, delimiter='\t')
	for row in reader:
		if row[0] != "genome":
			read = row[0]
			genome = map[row[0]]												# Foreach contig, map contig on genome
			contigs[row[0]][int(row[1])] = [int(row[2]), int(row[3])]			# new contig added, saved covered and total
cann = set()
for curr in contigs:
	genome = map[curr]
	for covered in contigs[curr]:	
		if covered in mappedcontigs[genome]:
			oldcov = mappedcontigs[genome][covered][0] 
			oldtot = mappedcontigs[genome][covered][1]
			mappedcontigs[genome][covered] = [oldcov + contigs[curr][covered][0], oldtot + contigs[curr][covered][1]]
		else:
			mappedcontigs[genome][covered] = [contigs[curr][covered][0], contigs[curr][covered][1]]
# scorro ogni contig, assegno un genoma e sommo i valori

with open("mapped_G18282_coverage.tsv","w") as tsv:
	for genome in mappedcontigs:
		for contig in mappedcontigs[genome]:
			cov=mappedcontigs[genome][contig][0]
			tot=mappedcontigs[genome][contig][1]
			freq=float(cov)/tot
			tsv.write('%s\t%i\t%i\t%i\t%f\n' % (genome, contig, cov, tot, freq))