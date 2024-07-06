# PSI-BLAST of consensus sequence against nr database with required parameters
# qcovs -> query coverage; pident -> percentage identity; length -> alignment length; slen -> target/reference seq length
import os
import sys
import time
import signal
from Bio import SeqIO
from Bio.Blast import NCBIWWW
from bs4 import BeautifulSoup as bs
import pandas as pd

def signal_handler(sig, frame):
    exit()
    
signal.signal(signal.SIGTERM, signal_handler)

def psiblast(dufID, pfamID, c):
    if int(c) == 0:
        inputFile = f"../DUFs/{dufID}/{pfamID}.fasta"
        outputFile = f"../DUFs/{dufID}/{pfamID}_psiblast.out"
    else: 
        inputFile = f"../DUFs/{dufID}/{pfamID}_Split_{c}.fasta"
        outputFile = f"../DUFs/{dufID}/{pfamID}_Split_{c}_psiblast.out"
        
    try:
        cmd_psiblast = f'psiblast -query {inputFile} -db /home/user/disk/proteinDB/nr -outfmt "6 qseqid saccver qcovs pident length slen salltitles mismatch gapopen qstart qend sstart send evalue bitscore" -out {outputFile} -out_pssm ../DUFs/{dufID}/{pfamID}.pssm -save_each_pssm -num_threads 8'
        
        os.system(cmd_psiblast)
        status = "psiblast completed"

    except:
        status = "Error during psiblast run"