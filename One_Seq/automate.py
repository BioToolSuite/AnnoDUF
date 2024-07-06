import os
import sys
import csv
import glob
import time
import signal
import pandas as pd
from Bio import SeqIO

# From here all imports are py file function
import split
import psiblast
import annotation
import blast_hit_title
import annotation_links

# DUF_File = pd.read_csv("Only_1_Seq_in_DUF_1.csv")
# ntProc = []

# dufID = sys.argv[1]
# pfamID = sys.argv[2]

def signal_handler(sig, frame):
    exit()
    
signal.signal(signal.SIGTERM, signal_handler)

def annoPipeline(dufID, pfamID, i):
    print('** PSI-BLAST in progress.... **')
    start_time = time.time()
    p = psiblast.psiblast(dufID, pfamID, i)
    end_time = time.time()
    execution_time = end_time - start_time
    print('** PSI-BLAST process completed in {:.2f} minutes {:.2f} seconds **'.format(execution_time // 60, execution_time % 60))
    
    print('** Storing blast hit titles and creating dictionary of occurrences of words **')
    b = blast_hit_title.psiblast_out_process(dufID, pfamID, i)

    # annotation
    if b.__eq__("Word dictionary created using retrieved blast hit titles") is True:
        print('** Filtering the Putative hits **')
        a = annotation.annotate(dufID, pfamID, i)

        if a.__eq__("Putative annotations obtained") is True:
            print("*************** Annotation is in process *****************")
            al = annotation_links.ann(dufID, pfamID, i)
        else:
            print("Unable to do putative annotation")
    else:
        print("Unable to do Filtering to putative hits")
    
df = pd.read_csv("../DUFs_with_1_Sequence.csv")
for i in range(len(df)):
    dufID = df["dufID"][i]
    pfamID = df["pfamID"][i]
    print(f"******************* {dufID} PROCESSING ****************")

    if len(SeqIO.read(f"../DUFs/{dufID}/{pfamID}.fasta", "fasta").seq) > 500:
        print('** Splitting is in progress.... **')
        splitSeq = split.splitter(dufID, pfamID)
        for i in range(1, int(splitSeq)):
            annoPipeline(dufID, pfamID, i)
        flst = glob.glob(f"../DUFs/{dufID}/{pfamID}_Split_*_putative_annotations.txt")
        df_list = []
        for filename in sorted(flst):
            df_list.append(pd.read_csv(filename))
        full_df = pd.concat(df_list)
        full_df.to_csv(f"../DUFs/{dufID}/{pfamID}_putative_annotations.txt", index=False)
    else:
        print("Sequence Length is not more than 500")
        annoPipeline(dufID, pfamID, 0)

