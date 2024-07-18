import pandas as pd
dufANDpfamID = pd.read_csv("dufPfamIDs.csv")
pfamLst = list(dufANDpfamID["pfamID"])
extract_lines = False
with open("Pfam-A.full", encoding='latin-1', errors='ignore') as file:
    for line in file:
        if line.startswith("#=GF AC") and line.replace('#=GF AC','').strip().split(".")[0] in pfamLst:
            accID = line.replace('#=GF AC','').strip().split(".")[0]
            extract_lines = True
            gf_cc_lines = []
            continue
        elif extract_lines and line.startswith("#=GF CC"):
            gf_cc_lines.append(line.strip())
        elif line.startswith("#=GF SQ"):
            if extract_lines:
                print(f"DUF ID: {accID}")
                if gf_cc_lines != []:
                    gf_cc_lines =  " ".join([line.replace("#=GF CC", "").lstrip().strip() for line in gf_cc_lines])
                    dufANDpfamID.loc[dufANDpfamID['pfamID'] == accID, 'description'] = gf_cc_lines
                else: 
                    dufANDpfamID.loc[dufANDpfamID['pfamID'] == accID, 'description'] = ""
            extract_lines = False
dufANDpfamID.to_csv("New_Test_List_of_DUFs.tsv", sep="\t", index=False)