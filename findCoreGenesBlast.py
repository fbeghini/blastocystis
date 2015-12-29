from Bio import Blast
from Bio.Blast import NCBIStandalone
import os, argparse

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument("--blastDB")
	parser.add_argument("--blastRef")
	parser.add_argument("--blastEXE", default='/usr/bin/blastn')

	args = parser.parse_args()

	blastDb = args.blastDB
	blastRef = args.blastRef
	blastExe = args.blastEXE

	blast_out, error_info = NCBIStandalone.blastall(blastExe, 'blastn', blastDb, blastRef)