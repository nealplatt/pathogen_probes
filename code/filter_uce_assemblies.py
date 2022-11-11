import os
import sys 
from Bio import SeqIO

in_file    = sys.argv[1]
out_dir   = sys.argv[2]
cov_cutoff = float(sys.argv[3])
len_cutoff = int(sys.argv[4])

out_file ="{}/{}".format(out_dir, os.path.basename(in_file))

if not os.path.exists(out_dir):
    os.makedirs(out_dir)

try:
    fasta_sequences = SeqIO.parse(open(in_file),'fasta')
    filtered_fasta_sequences = []

    for fasta in fasta_sequences:
        name, sequence = fasta.id, str(fasta.seq)
        cov=float(name.split("_")[-1])
        seq_len = int(name.split("_")[3])
        if ( (cov > cov_cutoff) and (seq_len > len_cutoff) ):
            filtered_fasta_sequences.append(fasta)

    if len(filtered_fasta_sequences)>1:
        SeqIO.write(filtered_fasta_sequences, out_file, "fasta")
except:
    print("{} does not contain any sequences".format(in_file))
    


