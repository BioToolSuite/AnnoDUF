import os
import collections
import re
import sys
import pandas as pd

def psiblast_out_process(dufID, pfamID, c):
    psiblast_output_file = f"../DUFs/{dufID}/{pfamID}_Split_{c}_psiblast.out"
    allTitlesFile = f"../DUFs/{dufID}/{pfamID}_Split_{c}_all_titles.txt"
    wordDictFile = f"../DUFs/{dufID}/{pfamID}_Split_{c}_word_dict.txt"
    delListFile = f"../DUFs/{dufID}/{pfamID}_Split_{c}_del_list.txt"
    if c == 0:
        psiblast_output_file = f"../DUFs/{dufID}/{pfamID}_psiblast.out"
        allTitlesFile = f"../DUFs/{dufID}/{pfamID}_all_titles.txt"
        wordDictFile = f"../DUFs/{dufID}/{pfamID}_word_dict.txt"
        delListFile = f"../DUFs/{dufID}/{pfamID}_del_list.txt"
    
    cnt1 = 0
    cnt2 = 0

    # deleting empty lines, if any in the psiblast output file
    cmd = f"sed -i '/^$/d' {psiblast_output_file}"
    os.system(cmd)

    # counting number of lines in psiblast output....to check if no hits obtained (empty file)
    blast_cnt = os.popen(f"wc -l < {psiblast_output_file}").read()

    if int(blast_cnt.rstrip()) > 0:
        lastLine = f"sed '${{/Search has CONVERGED!/d;}}' {psiblast_output_file} > temFile.txt && mv temFile.txt {psiblast_output_file}"
        os.system(lastLine)
        psibRst = pd.read_csv(psiblast_output_file, sep="\t", header=None)
        try:
            titles = open(allTitlesFile, 'w')           
            for i in range(0, len(psibRst)):
                if float(psibRst[2][i]) >= 80 and float(psibRst[3][i]) >= 20:
                    cnt1 += 1
                    target_cov = float(psibRst[4][i]) / float(psibRst[5][i]) * 100
                    if 80 <= target_cov <= 120:
                        cnt2 += 1
                        y = psibRst[6][i].split(sep="<>")
                        x = [item.rstrip().split(sep="[")[0].strip() for item in y]
                        titles.write('\n'.join(x)+'\n')
            titles.close()
        except Exception as e:
            # print(e)
            status = "Error while retrieving hits titles"

        if cnt1 == 0:
            status = "No significant hits found with query coverage greater than 80% and % identity greater than 20%"

        elif cnt1 != 0 and cnt2 == 0:
            status = "No hits fulfill target coverage 80-120%, although query coverage and percent identity criteria is fulfilled"

        else:
            # print("**Creating word dictionary**")
            file_path = allTitlesFile
            
            # Creating word dictionary
            word_dict = collections.OrderedDict()

            filter_word_list = []
            del_list = [dufID, pfamID]

            chars = ['(', ')', ':', '=', ',', "'", ';']

            # deleting unwanted characters from titles file
            for i in chars:
                cmd = f'sed -i "s/{i}//gI" {allTitlesFile}'
                os.system(cmd)

            # store count of words
            with open(file_path) as hits:
                for line in hits:
                    word_count(line.strip().split(' '), word_dict)

            word_list = list(word_dict.items())

            # list of common words to be eliminated during filtration for putative hits
            common_word_list = ['DUF', 'domain-containing', 'uncharacterized', 'uncharacterised', 'multispecies',
                                'wiith', 'to', 'of', 'in', 'similar', 'from', 'plasmid', '(plasmid)', 'contains',
                                'containing', 'partial', 'unknown', 'unnamed', 'hypothetical','RecName FullUncharacterized',
                                'RecName Full Uncharacterize', 'PREDICTED', 'LOW QUALITY PROTEIN']

            # filtering word list based eliminating low occurrences
            for i in word_list:
                occ = os.popen(f'grep -i -c "{i[0]}" {allTitlesFile}').read()
                if occ != '' and int(occ) > 5:
                    filter_word_list += [i]
                else:
                    del_list = del_list + [i[0]]

            # filtering word list by eliminating words in common word list
            for x in common_word_list:
                filter_word_list = [i for i in filter_word_list if not i[0].lower().startswith(x.lower())]
                del_list = del_list + [x]

            # filtering word list by eliminating words which are combination of letters and numbers

            filter_word_list = [i for i in filter_word_list if
                                not re.search(r'\d', i[0]) is not None and re.search(r'[a-zA-Z]', i[0]) is not None]
            del_list = del_list + [i[0] for i in filter_word_list if
                                   re.search(r'\d', i[0]) is not None and re.search(r'[a-zA-Z]', i[0]) is not None]

            # writing filtered word dict to file:
            with open(wordDictFile, "w") as f:
                for kv in filter_word_list:
                    f.write(str(kv) + '\n')
            f.close()

            # writing deleted words to file:
            with open(delListFile, "w") as f:
                for word in del_list:
                    f.write(str(word) + '\n')
            f.close()

            status = "Word dictionary created using retrieved blast hit titles"

            if len(filter_word_list) == 0:
                status = "Hits with qcov > 80% and pident > 20% are insignificant"
    else:
        status = "No hits obtained in psiblast run"

    return status

# Function to count total occurrences of word in entire file
def word_count(input_word, word_dict):
    for words in input_word:
        if words != '':
            if words.lower() in word_dict:
                word_dict[words.lower()] += 1
            else:
                word_dict[words.lower()] = 1
