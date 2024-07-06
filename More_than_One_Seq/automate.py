import os
import csv
import sys
import time
import pandas as pd

# From here all imports are py file function
import seqlen_filter
import msa
import consensus
import psiblast
import blast_hit_title
import annotation
import annotation_links

def main():
    fileName = sys.argv[1]
    
    df = pd.read_csv(f"../{fileName}")
    for i in range(len(df)):
        dufID = df["dufID"][i]
        pfamID = df["pfamID"][i]
        
        print(f"******************* {dufID} PROCESSING ****************")
    
        print('**Median seqlen filter**')
        s = seqlen_filter.seqlen_filter(f'../DUFs/{dufID}/{pfamID}.fasta')

        if s.__eq__("Median filter successful") is True:

            # Performing MSA of selected 40 sequences
            print('**Performing MSA**')
            m = msa.msa_clustal(f"../DUFs/{dufID}/{pfamID}_seqLen.fasta")

            # Generating consensus sequence from MSA
            print('** Generating consensus sequence from MSA **')
            c = consensus.cons(f"../DUFs/{dufID}/{pfamID}_Result.msa")

            # PSI-BLAST of consensus sequence against nr database
            print('** PSI-BLAST in progress.... **')
            start_time = time.time()
            p = psiblast.psiblast(f"../DUFs/{dufID}/{pfamID}_consensus.fasta")
            end_time = time.time()
            execution_time = end_time - start_time
            print('** PSI-BLAST process completed in {:.2f} minutes {:.2f} seconds **'.format(execution_time // 60, execution_time % 60))

            print('** Storing blast hit titles and creating dictionary of occurrences of words **')
            b = blast_hit_title.psiblast_out_process(f"../DUFs/{dufID}/{pfamID}_psiblast.out")

            # annotation
            if b.__eq__("Word dictionary created using retrieved blast hit titles") is True:
                print('** Filtering the Putative hits **')
                a = annotation.annotate(f"../DUFs/{dufID}/{pfamID}_all_titles.txt")

                if a.__eq__("Putative annotations obtained") is True:
                    print('** Annotating  the Putative hits **')
                    al = annotation_links.ann(f"../DUFs/{dufID}/{pfamID}_putative_hits.txt")

if __name__ == "__main__":
    main()
