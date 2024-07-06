import os
import sys
import random
import statistics
from Bio import SeqIO

# dufID = sys.argv[1]
# pfamID = sys.argv[2]
    
def seqlen_filter(fastaFile):
    dufANDpfam = fastaFile.split("/")
    dufID = dufANDpfam[2]
    pfamID = dufANDpfam[3].split(".")[0]
    try:
        length = []
        seq = []
        median_filter = []
        
        for record in SeqIO.parse(fastaFile, "fasta"):
            length+=[len(record.seq)]
            seq+=[record.format('fasta')]
            
        median = statistics.median(length)
        for i in range(len(length)):
            if median-20 <= length[i] <= median+20:
                median_filter.append(seq[i])

        if(len(median_filter))==0:
            status='Length of sequences not in significant range'

        elif len(median_filter) <= 40:
            with open(f"../DUFs/{dufID}/{pfamID}_seqLen.fasta", "w") as out:
                for i in median_filter:
                    out.write(i)
            status = "Median filter successful"
        
        else: 
            rand_nums = random.sample(range(0, len(median_filter)-1), 40)
            with open(f"../DUFs/{dufID}/{pfamID}_seqLen.fasta", "w") as out:
                for i in rand_nums:
                    out.write(median_filter[i])
            status = "Median filter successful"
    except Exception as e:
        status = "Error while filtering input sequences based on criteria"
        
    return status

# seqlen_filter(f"fileUpload/{dufID}/{pfamID}.fasta")