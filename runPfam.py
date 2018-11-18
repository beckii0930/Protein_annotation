import sys
import os

# get all the ac used for tool
def getData(input_dir):
	fh = open(input_dir,'r')
	Data = fh.read()
	Data = Data.split('>embl-cds:')[1:]
	Protein_ac = []
	Seq = ''

	# check out put dir for pfam raw output
	out_dir = "Pfam_out"
	if not os.path.exists(out_dir):
		os.makedirs(out_dir)

	for i in range (len(Data)):
		line = Data[i]
		newline = '>embl-cds:'
		ac = line.split()[0]
		# Seq=(newline + line)

		Protein_ac.append(ac)
		# f_out = open("Pfam_out/"+ ac + ".fasta", "w+")
		# f_out.write(Seq)
	return Protein_ac

def runPfam():
	if (len(sys.argv) < 1):
		print "Usage: python runPfam infile.fasta"

	# get sub sequences
	input_dir = sys.argv[1]
	Protein_ac = getData(input_dir)
	out_dir = "Pfam_out/"
	fasta_dir = "Fasta_files/"

	for ac in Protein_ac:
		# print ac
		command = "./hmmer-3.1b2-macosx-intel/binaries/hmmscan  -o "
		command += out_dir +ac + "_out.txt --noali Data/Pfam-A.hmm "
		command += fasta_dir +ac + ".fasta"
		os.system(command)



runPfam()