# Function to perform MSA of selected 40 sequences using clustal omega
import os
import sys

def msa_clustal(fastaSeqLen):
    dufANDpfam = fastaSeqLen.split("/")
    dufID = dufANDpfam[2]
    pfamID = dufANDpfam[3].replace('_seqLen', '').split(".")[0]
    try:
        cmd_msa = f"clustalo -i {fastaSeqLen} -o ../DUFs/{dufID}/{pfamID}_Result.msa --force"
        os.system(cmd_msa)
        status = 'MSA done'
    except:
        status = "Error while generating MSA"
        # print("Error while generating MSA")
    return status

# fastaFile = sys.argv[1]
# msa_clustal(fastaFile)
# msa_clustal("DUF6747/PF20532_seqLen.fasta")