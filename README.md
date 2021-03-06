# Protein_annotation
This is a protein annotation tool. It was first designed to annotate bacteria Acinetobacter baumannii in order to find out which proteins are potentially involved in mechanisms of antibiotic resistance.

The tool pipeline was written primarily in Python 3, and queries its data from databases such as BLAST, InterPro, Pfam, and Prosite. Each of these databases provides different information about a protein sequence. 

In order to make it more user-friendly, users are allowed to specify which tools to run, what their desired cutoff parameters are, and the output format (tab-, comma-, or space-separated plain text files). The tool pipeline can also take either FASTA-formatted sequence files or accession numbers as input, so that the user is not limited to only bacterial species.

## System Requirements:
To use this tool, it requires that you have all the following databases downloaded correctly. (See more details in project report)

## Usage:
python analysis.py [- flags]
-f infile.fasta
-a accession
-outTab
-outSpace
-outTsv
-noBlast
-noPfam
-noInterpro

## 1.	BLAST:
BLAST (Basic Local Alignment Search Tool) is a tool to infer functional and evolutionary relationships via alignment, as well as identify members of gene families. Since we were interested in finding possible A. baumannii protein functions, we used the NCBI BLASTp to search against protein sequence databases. We utilized the Biopython package in order to invoke the NCBI BLAST server over the internet, as well as use the built-in parsers in order to work with the raw XML files output by BLAST. We then wrote a Python script (blast.py) which queried our 100 protein sequences of interest and returned the output (which contained up to 10 hits) in two different file formats per sequence: a raw XML file output by BLAST, and a separate plain text file which contained the parsed information such as the protein accession number, gi number, predicted/known protein function, and E-value.
## 2.	InterPro:
As an integrated database, the InterPro consortium provides functional analysis of protein sequences by classifying them into families and predicting the presence of domains and other important sites. InterPro provides a software package known as InterProScan that allows sequences to be searched against key signatures in its database. Due to memory and processor limitations, we used InterPro’s REST service via their provided Python 3 client in order to query our sequences against their database over the internet. The drawback of this was that the InterPro REST service takes only one FASTA-formatted sequence at a time, and has a 30 sequence per batch limitation. We wrote a simple Python script (fasta_fragment.py) to fragment our file of 100 A. baumannii protein sequences into individual FASTA files. We also followed an online tutorial for the InterPro REST service to write a bash script to automate the input of up to 30 files at a time.1 Raw output was given as tab-separated plain text files, which we then parsed for protein function and relevant Gene Ontology terms.
## 3.	Pfam:
Pfam is a database of protein families and clans. It predicts the possible families the protein sequence belongs to using a HMM model. Since the web interface version of Pfam uses hmmscan to search, I downloaded the local version of HMMER 3.1b2. The database used by hmmscan tool is Pfam-A.hmm that can be downloaded from the ftp server but hmmpress should be run on the .hmm database. runPfam.py is used to query all sequences in the database. We used the default E-value of 1. It can be run independently using “python runPfam.py input_Dir”.
## 4.	Prosite:
Prosite is a database containing the function of proteins using regular expression. The perl script ps_scan.pl is downloaded and used to run Prosite search. The -d flag is used to specify the database file Prosite.dat that can be downloaded from the ftp server. The -s flag is used to mask all unspecific functions, such as N-glycosylation sites that are not crucial for studying the function of our protein sequences that can result in antibiotic resistance. -l flag is used to indicate level of sensitivity. 0 indicates only output results with high sensitivity while -1 indicates output results with low sensitivity. Since while running with high sensitivity doesn’t generate much output, I first query the database with high sensitivity (-l 0) and if there is no output, I run ps_scan again using low sensitivity (l=-1) to get more hits. runProsite.py is the code that runs Prosite to make it run, a gfortran compiler should be downloaded in the system. It can be run independently “python runProsite.py input_Dir”.
## 5.	Analysis:
Analysis.py is the python script that runs all tools together, concatenating and outputting all the results in one file. Users will run this script using the appropriate flags.
