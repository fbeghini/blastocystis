from Bio import Entrez

Entrez.email = "asfddsad@gmail.com"

ids = map(str.strip, open("/tmp/18SBlasto").readlines())
ids = [x.split("|")[1] for x in ids]

hpost = Entrez.epost("nuccore", id=','.join(ids))
postrecord = Entrez.read(hpost, validate=False)
query_key = postrecord["QueryKey"]
webenv = postrecord["WebEnv"]
hpost.close()
handle = Entrez.efetch(db="nucleotide", query_key=query_key, webenv=webenv)
result = Entrez.read(handle, validate=False)