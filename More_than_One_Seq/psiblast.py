# PSI-BLAST of consensus sequence against nr database with required parameters
# qcovs -> query coverage; pident -> percentage identity; length -> alignment length; slen -> target/reference seq length
import os
import sys

def psiblast(consensus_seq_file):
    dufANDpfam = consensus_seq_file.split("/")
    dufID = dufANDpfam[2]
    pfamID = dufANDpfam[3].replace('_consensus', '').split(".")[0]
    try:
        cmd_psiblast = f'psiblast -query {consensus_seq_file} -db /home/user/disk/proteinDB/nr -outfmt "6 qseqid saccver qcovs pident length slen salltitles mismatch gapopen qstart qend sstart send evalue bitscore" -out ../DUFs/{dufID}/{pfamID}_psiblast.out -out_pssm ../DUFs/{dufID}/{pfamID}.pssm -save_each_pssm -num_threads 8'
        
        os.system(cmd_psiblast)
        status = "psiblast completed"

    except:
        status = "Error during psiblast run"

    return status