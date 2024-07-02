import pandas as pd
df = pd.DataFrame(columns=["dufID", "pfamID"])
c = 0
lastLine = ""
with open("Pfam-A.full", encoding='latin-1', errors='ignore') as fileobject:
    for line in fileobject:
        if line.startswith('#=GF ID') and "DUF" in line:
            ID = line.replace('#=GF ID','').strip()
            c = c + 1
        elif line.startswith('#=GF AC') and c == 1:
            acc = line.replace('#=GF AC','').strip()
            df.loc[len(df)] = ID, acc
            c = 0
            print(f"DUF ID: {ID}, Pfam Accession ID: {acc}")
        else:
            lastLine = line
df.to_csv("List_of_DUFs_along_PfamID.tsv", sep="\t", index=False)