import sys
from Bio import SeqIO

in_fasta = sys.argv[1]

good_seqs = [] 
with open(in_fasta, "r") as in_file: 
     for record in SeqIO.parse(in_file, "fasta"): 
         if (len(record.seq) == 120) and (record.seq.count("N") == 0): 
             good_seqs.append(record) 

with open("repeat_cleaned_seqs.fas", "w") as output_file:
    SeqIO.write(good_seqs, output_file, "fasta")



