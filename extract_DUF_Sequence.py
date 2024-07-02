import os
import pandas as pd
from Bio import SeqIO
from multiprocessing import Pool, Manager

def process_row(queue):
    while not queue.empty():
        row = queue.get()
        if row is None:
            break
        dufID = row["dufID"]
        pfamID = row["pfamID"]
        os.system(f"bash extract_matching_fasta.sh Pfam-A.fasta {pfamID} {dufID}")
        print(f"DUF ID: {dufID}")

if __name__ == "__main__":
    dufpfamFile = pd.read_csv("List_of_DUFs_along_PfamID.tsv", sep="\t")
    
    # Create a Manager Queue
    manager = Manager()
    queue = manager.Queue()

    # Put rows into the queue
    for index, row in dufpfamFile.iterrows():
        queue.put(row)

    # Create a pool of workers with 20 processes
    with Pool(20) as pool:
        # Use the pool to apply the process_row function to the queue
        pool.map(process_row, [queue for _ in range(20)])