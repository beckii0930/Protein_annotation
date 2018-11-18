# import runPfam
# import runProsite
import sys
import os
import csv
import argparse

def main():
	if (len(sys.argv) < 1):
		print "Usage: python analysis.py"
		print "-f infile.fasta"
		print "-a accession"
		print "-outTab"
		print "-outSpace"
		print "-outTsv"
		print "-noBlast"
		print "-noPfam"
		print "-noInterpro"

	# # create single fasta files from seqeunce fasta files
	# input_dir = sys.argv[1]

	# run Pfam adn prosite on the fasta files
	#python analysis.py  subSeq.fasta

	parser = argparse.ArgumentParser()
	parser.add_argument('-f', dest='input_dir')
	parser.add_argument('-a', dest='ac')
	parser.add_argument('-outTab', action='store_true', dest='outTab')
	parser.add_argument('-outSpace', action='store_true', dest='outSpace')
	parser.add_argument('-noBlast', action='store_true', dest='noBlast')
	parser.add_argument('-noPfam', action='store_true', dest='noPfam')
	parser.add_argument('-noProsite', action='store_true', dest='noProsite')
	args = parser.parse_args()
	# print args
	if args.ac is not None:
		print "value of accession is " + str(args.ac)
		if args.noPfam is not None:
			if args.noPfam == False:
				pfam_command =  "python runPfam.py --acc " + str(args.ac)
				print "----- running pfam"
				# os.system(pfam_command)
		if args.noProsite is not None:
			if args.noProsite == False:
				ps_command =  "python runProsite.py -f " + str(args.ac)
				print "----- running prosite"
				# os.system(ps_command)
		Protein_ac = [args.ac]

	elif args.input_dir is not None:
		Protein_ac = getData(args.input_dir)
		input_dir = args.input_dir
		if args.noBlast is not None:
			if args.noBlast == False:
				blast_command = "python blast.py -f " + input_dir
				print "----- running BLAST"
				# os.system(blast_command)

		if args.noPfam is not None:
			if args.noPfam == False:
				pfam_command =  "python runPfam.py " + input_dir
				print "----- running pfam"
				# os.system(pfam_command)

		if args.noProsite is not None:
			if args.noProsite == False:
				pfam_command =  "python runProsite.py " + input_dir
				print "----- running prosite"
				# os.system(ps_command)
	else:
		print 'Error: no input specified!'
		return;

	# integrate data
	delimiter = ''
	if args.outSpace == True:
		delimiter = ' '
	elif args.outTab == True:
		delimiter = '\t'

	# integrateData(Protein_ac, delimiter)
	# integrateInterpro(Protein_ac, delimiter)

def integrateInterpro(Protein_ac, delimiter):
	out_dir = "Results"
	# out_dir = "interpro_out/"

	if not os.path.exists(out_dir):
		os.makedirs(out_dir)
	f_out = open(out_dir + "/output_interpro.txt", 'w+')
	header = 'AC' +delimiter + 'Prosite' + delimiter + 'Pfam' + delimiter + 'GO\n'
	line = header
	pfam_dict, prosite_dict, go_dict = interpro_other(Protein_ac)

	# print pfam_dict
	for ac in Protein_ac:
		line +=ac + delimiter
		interpro_start = True
		if (ac not in pfam_dict.keys()):
			line += "No pfam Hit"
		else:
			for m in range (len(pfam_dict[ac])):
				# for n in range (len(pfam_dict[ac][m])):
				if interpro_start:
					interpro_start = False
					line += pfam_dict[ac][m][0] + '/' + pfam_dict[ac][m][1]

				else:
					line += "," + pfam_dict[ac][m][0] + '/' + pfam_dict[ac][m][1]

		line += delimiter

		interpro_start_2 = True
		if (ac not in prosite_dict.keys()):
			line += "No prosite Hit"
		else:
			for n in range (len(prosite_dict[ac])):
				if interpro_start_2:
					interpro_start_2 = False
					line += prosite_dict[ac][n][0] + '/' +  prosite_dict[ac][n][1]
				else:
					line += "," + prosite_dict[ac][n][0] + '/' +  prosite_dict[ac][n][1]

		line += delimiter
		interpro_start_3 = True
		if (ac not in go_dict.keys()):
			line += "No GO Hit"
		else:
			for i in range (len(go_dict[ac])):
				if interpro_start_3:
					interpro_start_3 = False
					line += go_dict[ac][i]

				else:
					line += "," + go_dict[ac][i]
		line += '\n'
	f_out.write(line)



def integrateData(Protein_ac, delimiter):
	out_dir = "Results"
	if not os.path.exists(out_dir):
		os.makedirs(out_dir)
	f_out_others = open(out_dir + "/output_others.txt", 'w+')
	f_out_ours = open(out_dir + "/output_ours.txt", 'w+')
	header = 'AC' + delimiter + 'Pfam' + delimiter + 'Prosite' + delimiter + 'BLAST' + delimiter + 'Interpro\n'
	line_others = header
	line_ours = header
	Prosite_out_dir = "Prosite_out/"
	Pfam_out_dir = "Pfam_out/"

	#get the prosite keywords
	Prosite_out_data,Prosite_out_data_f =  getPrositeKey(Protein_ac)

	#get the pfam keywords
	Pfam_out_data = getPfamKey(Protein_ac)
	blast_out_data = blast_to_dict(Protein_ac)
	interpro_out_data = interpro_parse(Protein_ac)
	# f_out.write(header)

	print interpro_out_data['ABO10621']

	for ac in Protein_ac:
		line_others += ac + delimiter
		line_ours += ac + delimiter
		ps_start = True
		pfam_start =  True
		blast_start = True
		interpro_start = True

		##----------- prosite##----------- 
		for i in range (len(Prosite_out_data[ac])):
			if ps_start:
				ps_start = False
				line_others += Prosite_out_data[ac][i][0]+"/"+ Prosite_out_data[ac][i][2]
				# line_ours += Prosite_out_data[ac][i][2]

			else:
				line_others += "," + Prosite_out_data[ac][i][0]+"/"+ Prosite_out_data[ac][i][2]
				# line_ours += "," + Prosite_out_data[ac][i][2]

		if (len(Prosite_out_data[ac]) == 0):
			line_others += "No Prosite Hit"
			# line_ours += "No Prosite Hit"
		line_others += delimiter
	 	# line_ours += '\t'

	 	ps_start = True
	 	for i in range (len(Prosite_out_data_f[ac])):
			if ps_start:
				ps_start = False
				# line_others += Prosite_out_data[ac][i][0]+"/"+ Prosite_out_data[ac][i][2]
				line_ours += Prosite_out_data_f[ac][i][2]

			else:
				# line_others += "," + Prosite_out_data_f[ac][i][0]+"/"+ Prosite_out_data_f[ac][i][2]
				line_ours += "," + Prosite_out_data_f[ac][i][2]

		if (len(Prosite_out_data_f[ac]) == 0):
			# line_others += "No Prosite Hit"
			line_ours += "No Prosite Hit"
		# line_others += '\t'
	 	line_ours += delimiter

		##----------- pfam##----------- 
		for j in range (len(Pfam_out_data[ac])):
			if pfam_start:
				pfam_start = False
				line_others +=Pfam_out_data[ac][j][0]+"/"+Pfam_out_data[ac][j][1]
				line_ours +=Pfam_out_data[ac][j][1]

			else:
				line_others += "," + Pfam_out_data[ac][j][0]+"/"+Pfam_out_data[ac][j][1]
				line_ours += "," + Pfam_out_data[ac][j][1]

	 	if (len(Pfam_out_data[ac])==0):
			line_others += "No Pfam Hit"
			line_ours += "No Pfam Hit"

		line_others += delimiter
	 	line_ours += delimiter


	 	##----------- blast##----------- 
		for k in range (len(blast_out_data[ac])):
			if blast_start:
				blast_start = False
				line_others += blast_out_data[ac][k][0]+"/"+ blast_out_data[ac][k][1]+"/"+ blast_out_data[ac][k][2]
				line_ours += blast_out_data[ac][k][2]

			else:
				line_others += "," + blast_out_data[ac][k][0]+"/"+ blast_out_data[ac][k][1]+"/"+ blast_out_data[ac][k][2]
				line_ours += "," + blast_out_data[ac][k][2]

		if (len(blast_out_data[ac]) == 0):
			line_others += "No blast Hit"
			line_ours += "No blast Hit"
	 	line_others += delimiter
	 	line_ours += delimiter

	 	##----------- interpro##----------- 
	 	if (ac not in interpro_out_data.keys()):
			line_others += "No interpro Hit"
			line_ours += "No interpro Hit"
		else:
			for m in range (len(interpro_out_data[ac])):
				if interpro_start:
					interpro_start = False
					line_others += interpro_out_data[ac][m]
					line_ours += interpro_out_data[ac][m]

				else:
					line_others += "," + interpro_out_data[ac][m]
					line_ours += "," + interpro_out_data[ac][m]


	 	line_others += '\n'
	 	line_ours += '\n'
	f_out_others.write(line_others)
	f_out_ours.write(line_ours)


def getPrositeKey(Protein_ac):
	print "in get prosite key"
	# Protein_ac = ['ABO10513', 'ABO10514']
	Data = {}
	Data_f={}
	Prosite_out_dir = "Prosite_out/"
	for ac in Protein_ac:
		f_in = open(Prosite_out_dir + ac + "_out.txt", 'r')
		f_in_f = open(Prosite_out_dir + ac + "_out_filtered.txt", 'r')
		raw_data = f_in.read().split('>embl-cds:')
		raw_data_f = f_in_f.read().split('>embl-cds:')
		Data[ac] = []
		Data_f[ac] = []
		# get the id for the tools
		for i in raw_data:
			if ac in i:
				# print i
				hit = (i.split(".")[0]).split()
				id = hit[2]
				key = hit[3]
				description = ' '.join(hit[4:])
				# print ac , id, key, description
				Data[ac].append([id, key, description])
		for i in raw_data_f:
			if ac in i:
				hit = (i.split(".")[0]).split()
				id = hit[2]
				key = hit[3]
				description = ' '.join(hit[4:])
				# print ac , id, key, description
				Data_f[ac].append([id, key, description])
	# print Data
	return Data, Data_f

def getPfamKey(Protein_ac):
	print "in get prosite key"
	# Protein_ac = ['ABO10513', 'ABO10514']
	Data = {}
	Pfam_out_dir = "Pfam_out/"
	for ac in Protein_ac:
		f_in = open(Pfam_out_dir + ac + "_out.txt", 'r')
		raw_data = f_in.readlines()
		Data[ac] = []
		# get the id for the tools
		for i in raw_data:
			if ">>" in i:
				# print i
				hit = i.split()
				key = hit[1]
				description = ' '.join(hit[2:])
				# print ac , key, description
				Data[ac].append([key, description])
	return Data

def blast_to_dict(Protein_ac):
	blast_dict = {}
	blast_file = "blast_out/"
	for ac in Protein_ac:
		with open(blast_file+ac+ "_out.txt", 'r') as raw_blast:
			# print('\n')
			# print('FILENAME:', blast_file)
			for line in raw_blast:
				if line.startswith('gi:'):
					gi = line.strip()
					gi = gi.split(':')[1]
					accession = next(raw_blast).strip()
					accession = accession.split(':')[1]

					function = next(raw_blast).strip()
					if 'MULTISPECIES' in function:
						split = function.split(':')
						function = split[1] + ':' + split[2]
					else:
						function = function.split(':')[1]

					if ac in blast_dict.keys():
						blast_dict[ac].append([accession, gi, function])
					else:
						blast_dict[ac] = []
						blast_dict[ac].append([accession, gi, function])
	return blast_dict

def interpro_parse(Protein_ac):
	file_path = 'interpro_out/'
	# ac = interpro_file.strip('.tsv.txt')
	go_dict = {}
	for ac in Protein_ac:
		with open(file_path + ac + '.tsv.txt', 'r') as infile:
			tsv = csv.reader(infile, delimiter = '\t')
			for row in tsv:
				# Extract GO terms
				if len(row) > 13:
					if row[13].startswith('GO'):
					#print(row[0] + '\t' + row[13])

						go_list = row[13].split('|')
						#print(go_list)

						if ac in go_dict.keys():
							for go in go_list:
								if go not in go_dict[ac]:
									go_dict[ac].append(go)
						else:
							go_dict[ac] = []
							for go in go_list:
								if go not in go_dict[ac]:
									go_dict[ac].append(go)
			# if len(go_dict) > 0:
				# print(go_dict)
	return go_dict


def interpro_other(Protein_ac):
	file_path = 'interpro_out/'
	# ac = interpro_file.strip('.tsv.txt')
	go_dict = {}
	pfam_dict = {}
	prosite_dict = {}
	for ac in Protein_ac:
		with open(file_path + ac + '.tsv.txt', 'r') as infile:
			tsv = csv.reader(infile, delimiter = '\t')
			for row in tsv:
				# Extract Pfam ID and function description
				if 'Pfam' in row:
					pfam_id = row[4]
					pfam_desc = row[5]

					if ac in pfam_dict.keys():
						pfam_dict[ac].append([pfam_id, pfam_desc])
					else:
						pfam_dict[ac] = [[pfam_id, pfam_desc]]

				# Extract ProSiteProfiles ID and function description
				if 'ProSiteProfiles' in row:
					prosite_id = row[4]
					prosite_desc = row[5]

					if ac in prosite_dict.keys():
						prosite_dict[ac].append([prosite_id, prosite_desc])
					else:
						prosite_dict[ac] = [[prosite_id, prosite_desc]]

				# Extract GO terms
				if len(row) > 13:
					if row[13].startswith('GO'):
					#print(row[0] + '\t' + row[13])

						go_list = row[13].split('|')
						#print(go_list)

						if ac in go_dict.keys():
							for go in go_list:
								if go not in go_dict[ac]:
									go_dict[ac].append(go)
						else:
							go_dict[ac] = []
							for go in go_list:
								if go not in go_dict[ac]:
									go_dict[ac].append(go)
	return pfam_dict, prosite_dict, go_dict

def getData(input_dir):
	fh = open(input_dir,'r')
	Data = fh.read()
	Data = Data.split('>embl-cds:')[1:]
	Seq = ''
	Protein_ac =[]
	# check out put dir for pfam raw output
	out_dir = "Fasta_files"

	if not os.path.exists(out_dir):
		os.makedirs(out_dir)

	for i in range (len(Data)):
		line = Data[i]
		newline = '>embl-cds:'
		ac = line.split()[0]
		Seq=(newline + line)

		Protein_ac.append(ac)

		f_out = open("Fasta_files/"+ ac + ".fasta", "w+")
		f_out.write(Seq)
	return Protein_ac

if __name__ == "__main__":
    main()