import random 
import sys
from Bio import SeqIO 
from Bio.Seq import Seq
 
num_passes = int(sys.argv[1]) 
in_fas     = sys.argv[2] 
out_fas    = sys.argv[3]
 
seq_records = list(SeqIO.parse(in_fas, "fasta")) 
 
randomized_seqs=[]  
 
with open(out_fas, 'w') as out_file:
    for seq in seq_records: 
        #make mulitple passes to randomize nucleotide order 
        for i in list(range(1,num_passes+1)): 
            seq.seq = Seq(''.join(random.sample(str(seq.seq),len(seq.seq)))) 
        
        #write to file
        seq = ">{}\n{}\n".format(seq.seq, seq.description)
        out_file.write(seq)
