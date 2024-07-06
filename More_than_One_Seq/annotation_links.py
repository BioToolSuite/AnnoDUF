import os
import sys
import pandas as pd
import requests
from bs4 import BeautifulSoup

proxies = {'https':'http://proxy.ibab.ac.in:3128', 'http':'http://proxy.ibab.ac.in:3128'}
def uniprot(searchText):
    searchText = searchText.replace(" ", "%20")
    url = f'https://rest.uniprot.org/uniprotkb/search?query={searchText}'
    data = requests.get(url, proxies=proxies).json()
    if "messages" in data or len(data['results']) == 0:
        return "No hits!"
    else:
        return f"https://www.uniprot.org/uniprotkb?query={searchText}"
    
def Convert(string):
    li = list(string.split("\n"))
    return li

def kegg_Pathway(searchText):
    form_data = {}
    url = "https://www.genome.jp/kegg/pathway.html"
    proxies = {'https':'http://proxy.ibab.ac.in:3128/', 'http':'http://proxy.ibab.ac.in:3128/'}
    response = requests.get(url, proxies=proxies)
    soup = BeautifulSoup(response.content, 'html.parser')
    form = soup.find('form')
    form_data['keyword'] = searchText
    response = requests.post(form['action'], data=form_data, proxies=proxies)
    soup = BeautifulSoup(response.content, "html.parser")
    if Convert(str(soup))[0] == "Bad request":
        return "No hits!"
    else:
        return soup.find("form").select("table tr td font")[0].text.strip()
    
def ann(putHitFile):
    # print("*************** Annotation is in process *****************")
    dufANDpfam = putHitFile.split("/")
    dufID = dufANDpfam[2]
    pfamID = dufANDpfam[3].replace('_putative_hits', '').split(".")[0]
    with open(putHitFile, 'r') as h:
        hits = h.read().splitlines()
    h.close()
    
    # Remove unnecessary keyword from the list of annotations
    y = [i.strip() for i in sorted(hits)]
    putative_annotations = []
    b = [putative_annotations.append(i) for i in y if i.rstrip() not in putative_annotations]
    del_ann = ['predicted protein', 'putative protein', 'conserved protein', 'putative conserved protein', 'predicted conserved protein']
    del_contain_ann = ['RecName FullUncharacterized protein', 'RecName Full Uncharacterized protein', dufID, pfamID]
    def rmvUnwntKeyWD(delAnnp, putaAnno):
        commSet = list(set(delAnnp) & set(putaAnno))
        putaAnno = list(set(putaAnno) - set(commSet))
        return putaAnno
    putative_annotations = rmvUnwntKeyWD(del_ann, putative_annotations)
    putative_annotations = sorted(rmvUnwntKeyWD(del_contain_ann, putative_annotations))
    
    # sort the file and remove the duplicates doesn't matter it has special character or not
    my_list = []
    previous_value = None
    c = 0
    for pf in sorted(set(putative_annotations)):
        if previous_value is not None:
            if previous_value != pf:
                my_list.append(pf)
        previous_value = pf
        if c == 0:
            first_value = pf
        if len(putative_annotations)-1 == c and first_value not in my_list:
            my_list.append(first_value)
        c+=1
    my_list = sorted(set(my_list))
    
    # Create the final list of annotation
    df = pd.DataFrame(columns=["Putative_Annotation", "Uniprot", "Kegg"])
    for i in my_list:
        blkLst = []
        blkLst.append(i)
        blkLst.append(uniprot(i))
        if kegg_Pathway(i) != "No hits!":
            j = i.replace(" ", "+")
            blkLst.append(f"https://www.kegg.jp/kegg-bin/search_pathway_text?map=map&keyword={j}&mode=1&viewImage=true")
        else:
            blkLst.append("No hits")
        df.loc[len(df)] = blkLst
    df.to_csv(f'../DUFs/{dufID}/{pfamID}_putative_annotations.txt', index=False)
    # print("*************** Annotation is completed ******************")
    
# putHitFile = sys.argv[1]
# ann(putHitFile)

# ann("DUF1079/PF06435_putative_hits.txt")
# ann(putHitFile)
