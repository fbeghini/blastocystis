from Bio import Entrez
Entrez.email = "francesco.beghini@studenti.unitn.it"
handle = Entrez.esearch(db="nucleotide", term="Giardia", retmax="19684")
record = Entrez.read(handle)

for item in record["IdList"]:
	hitem = Entrez.efetch(db="nucleotide", id=item, rettype="fasta")
	print hitem.read()