""" 
Pipeline for running BLAST on a list of FASTA-formatted sequences
Usage:

Note: Due to NCBI server limitations, you should not run more than 50 sequences
at a time. It it best to split up your queries into multiple files if you have
many sequences to run. For more information please see [NCBI USAGE GUIDELINES].
"""

import os
import argparse
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import SeqIO

E_VALUE_THRESH = 0.04    # E-value threshold for alignments

# TODO Re-write this script to allow multiprocessing
# TODO write usage instructions; probably use argparse to accept either 
# FASTA or accessions

#------------------------------------------------------------------------------
# Argument parsing
parser = argparse.ArgumentParser()
parser.add_argument('filename', help='Path to the input file')
parser.add_argument('-f', help='Option to take FASTA file as input', 
                    action='store_true')
parser.add_argument('-a', help='Option to take list of accession numbers as input',
                    action='store_true')
args = parser.parse_args()
#------------------------------------------------------------------------------
# If the '-f' option is specified, take FASTA-formatted file as input
if args.f:
    # For each FASTA-formatted sequence read in from file as SeqRecord object
    for seq_record in SeqIO.parse(args.filename, 'fasta'):

        # Print the sequence description
        print("Sequence Description:", seq_record.description)
        print("Running BLAST...")

        # Run BLAST over internet and save top 10 hits
        result_handle = NCBIWWW.qblast('blastp', 'nr', 
                                       seq_record.format('fasta'), hitlist_size=10)

        print("BLAST done!")

        # Save a local copy of the BLAST output file
        protein_name = seq_record.description.split()[0].split(':')[1]
        print(protein_name)
        with open(os.path.join('C:/Users/Alyssa/Desktop/CSE182/project/blast_raw',
                  str(protein_name) + '_raw.xml'), 'w') as out_handle:
            out_handle.write(result_handle.read())
        result_handle.close()

        # Load the BLAST result back onto handle
        result_handle = open(os.path.join('C:/Users/Alyssa/Desktop/CSE182/project/blast_raw',
                  str(protein_name) + '_raw.xml'))

        # Parse the BLAST output
        blast_record = NCBIXML.read(result_handle)

        with open(os.path.join('C:/Users/Alyssa/Desktop/CSE182/project/blast_out', 
                  str(protein_name) + '_out.txt'), 'w') as outfile:
            for alignment in blast_record.alignments:
                # Can have more than one result due to protein redundancy        
                title = alignment.title
                gi = title.split('|')[1]
                accession = title.split('|')[3]
                function = title.split('|')[4]
                function = function.replace('>gi', '')

                outfile.write('gi:')
                outfile.write(gi)
                outfile.write('\n')

                outfile.write('accession:')
                outfile.write(accession)
                outfile.write('\n')

                outfile.write('function:')
                outfile.write(function)
                outfile.write('\n')
                outfile.write('\n')

#------------------------------------------------------------------------------
# If the '-a' option is specified, take list of accession numbers as input
# TODO after presentation

if args.a:
    pass