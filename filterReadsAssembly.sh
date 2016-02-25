#!/bin/bash
mkdir assembly_contigs
for x in ST4*.fna; do
	echo $x
	subtype=${x%_*}
	if [[ "$subtype" == "ST4" ]]
	then bdb="blastdb/GCF_000743755"
	else bdb=`ls blastdb/$subtype* | grep JZR | grep nin`
	fi
	blastn -outfmt 6 -query $x -db ${bdb%.*} -evalue 100 | cut -f1 | sort | uniq > assembly_contigs/$x.lst
	python ~/pyphlan/fna_sss.py $x --select --ids assembly_contigs/$x.lst > ${x%.*}.fna.final
done
