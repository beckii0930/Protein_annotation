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
	out_dir = "Prosite_out"
	if not os.path.exists(out_dir):
		os.makedirs(out_dir)

	for i in range (len(Data)):
		line = Data[i]
		newline = '>embl-cds:'
		ac = line.split()[0]
		# Seq=(newline + line)

		Protein_ac.append(ac)
		# f_out = open("Prosite_out/"+ ac + ".fasta", "w+")
		# f_out.write(Seq)
	return Protein_ac

def runProsite():
	if (len(sys.argv) < 1):
		print "Usage: python runProsite infile.fasta"

	# get sub sequences
	input_dir = sys.argv[1]
	Protein_ac = getData(input_dir)
	out_dir = "Prosite_out/"
	fasta_dir = "Fasta_files/"

	# Protein_ac = ['ABO10513']
	for ac in Protein_ac:
		print ac
		command = "perl ps_scan/ps_scan.pl -d Data/prosite.dat "
		command += fasta_dir + ac + ".fasta > " 
		command += out_dir + ac + "_out.txt"
		os.system(command)

		command = "perl ps_scan/ps_scan.pl -s -d Data/prosite.dat "
		command += fasta_dir + ac + ".fasta > " 
		command += out_dir + ac + "_out_filtered.txt"
		os.system(command)

		# if no hits fount, use a lower cutoff
		if os.stat(out_dir + ac + "_out_filtered.txt").st_size == 0:
			print ac + " has no high scoring hits, using a lower cut off"
			command = "perl ps_scan/ps_scan.pl -s -l -1 -d Data/prosite.dat "
			command += fasta_dir + ac + ".fasta > " 
			command += out_dir + ac + "_out_filtered.txt"
			os.system(command)
		if os.stat(out_dir + ac + "_out_filtered.txt").st_size == 0:
			print "Lowering the hits still no output"



runProsite()