"""
Script which, given a raw BLAST output file (.XML format), extracts the E-value data.
"""

import os
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import SeqIO

file_path = 'C:/Users/Alyssa/Desktop/CSE182/project/blast_raw/'

def get_evalues():
    with open('C:/Users/Alyssa/Desktop/CSE182/project/blast_eval.txt', 'w') as outfile:
        # Parse all raw BLAST output files in directory
        for rawblast_file in os.listdir(file_path):

            # Load the BLAST result back onto handle
            result_handle = open(file_path + rawblast_file)

            # Parse the BLAST output
            blast_record = NCBIXML.read(result_handle)
    
            for alignment in blast_record.alignments:
                # Can have more than one result due to protein redundancy
                
                for hsp in alignment.hsps:
                    e_value = hsp.expect
                    outfile.write(str(e_value) + ',')
            outfile.write('\n')

    # Strip the trailing comma in each row
    with open('C:/Users/Alyssa/Desktop/CSE182/project/blast_eval.txt', 'r') as readin:
        with open('C:/Users/Alyssa/Desktop/CSE182/project/FINAL_blast_eval.txt', 'w') as eval_out:
            for line in readin:
                line = line.rstrip(',\n')
                eval_out.write(str(line))
                eval_out.write('\n')

get_evalues()