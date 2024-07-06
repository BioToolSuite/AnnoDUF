# Generating consensus sequence from MSA
import os
import sys

def cons(msa_file):
    dufANDpfam = msa_file.split("/")
    dufID = dufANDpfam[2]
    pfamID = dufANDpfam[3].replace('_Result', '').split(".")[0]
    try:
        cmd_cons = f"cons -sequence {msa_file} -outseq ../DUFs/{dufID}/{pfamID}_consensus.fasta -name {dufID}_consensus"
        os.system(cmd_cons)
        status = "consensus sequence generated"
    except:
        status = "Error while generating consensus sequence from MSA"

    return status