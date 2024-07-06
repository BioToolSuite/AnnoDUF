from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def writeSeq(seq, frag_id, desc, locate, c):
    my_seqs = SeqRecord(Seq(seq), id=frag_id, description=desc)
    SeqIO.write(my_seqs, f"{locate}_Split_{c}.fasta", "fasta")
    
def splitter(dufID, pfamID):
    locate = f"../DUFs/{dufID}/{pfamID}"
    input_file = f"{locate}.fasta"
    c = 1
    aligned_fragments = []
    previousValue = None
    
    with open(input_file) as input_handle:
        records = list(SeqIO.parse(input_handle, "fasta"))
        for record in records:
            recordID = record.id
            recordDesc = record.description
            recordSeq = record.seq
            for i in range(0, len(recordSeq), 250):
                fragment_seq = recordSeq[i:i+500]
                fragment_id = f"{recordID}_{i+1}-{i+len(fragment_seq)}"
                if previousValue is not None:
                    if previousValue != i+len(fragment_seq):
                        writeSeq(fragment_seq, fragment_id, recordDesc, locate, c)
                        c += 1
                else:
                    writeSeq(fragment_seq, fragment_id, recordDesc, locate, c)
                    c += 1
                previousValue = i+len(fragment_seq)
    return c