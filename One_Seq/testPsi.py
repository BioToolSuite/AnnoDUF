import os
import sys
from Bio import SeqIO
from Bio.Blast import NCBIWWW
from bs4 import BeautifulSoup as bs
import pandas as pd


def psiblast(dufID, pfamID, c):
    print(dufID, pfamID, c)
    if int(c) == 0:
        inputFile = f"../../fileUpload/{dufID}/{pfamID}.fasta"
        outputFile = f"../../fileUpload/{dufID}/{pfamID}_psiblast.out"
    else: 
        inputFile = f"../../fileUpload/{dufID}/{pfamID}_Split_{c}.fasta"
        outputFile = f"../../fileUpload/{dufID}/{pfamID}_Split_{c}_psiblast.out"
        
    # dufANDpfam = consensus_seq_file.split("/")
    # dufID = dufANDpfam[1]
    # pfamID = dufANDpfam[2].replace('_consensus', '').split(".")[0]
    seq_record = next(SeqIO.parse(open(inputFile),'fasta')) 
    result_handle = NCBIWWW.qblast("blastp", "nr", seq_record.seq, format_type="HTML", hitlist_size=500, service ="psi")
    result = result_handle.read()
    soup = bs(result, 'html.parser')
    table = soup.find("table", {"id": "dscTable"}).find("tbody").findAll("tr")
    df = pd.DataFrame(columns=["desc", "queryCover", "percentIden", "seqLen", "acclen"])
    for row in table:
        desc = row.findAll("td")[1].find("a").text
        slen = row.findAll("td")[1].find("a")["len"]
        qCov = row.findAll("td")[7].text.rstrip("%")
        piden = row.findAll("td")[9].text.rstrip("%")
        acclen = row.findAll("td")[10].text
        df.loc[len(df.index)] = [desc, qCov, piden, slen, acclen]
    df.to_csv(outputFile, index=False)
    
import time
print('** PSI-BLAST in progress.... **')
start_time = time.time()
psiblast(751, "Seq_1", 0)
end_time = time.time()
execution_time = end_time - start_time
print('** PSI-BLAST process completed in {:.2f} minutes {:.2f} seconds **'.format(execution_time // 60, execution_time % 60))