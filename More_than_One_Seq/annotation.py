import os
import sys

def annotate(blastHitRstFile):
    dufANDpfam = blastHitRstFile.split("/")
    dufID = dufANDpfam[2]
    pfamID = dufANDpfam[3].replace('_all_titles', '').split(".")[0]
    
    # print("**Annotation under process**")
    # creating a temporary copy of hit->titles file
    cmd = f'sort -f {blastHitRstFile} | uniq -i > ../DUFs/{dufID}/{pfamID}_all_titles_tmp.txt'
    os.system(cmd)

    hit_words = []
    del_words = []

    # appending words from filtered word dictionary
    with open(f"../DUFs/{dufID}/{pfamID}_word_dict.txt") as f:
        for word in f:
            hit_words.append((word.rstrip().split("'"))[1])
    f.close()

    # appending words which were deleted to generate filtered word dictionary
    with open(f"../DUFs/{dufID}/{pfamID}_del_list.txt") as f:
        for word in f:
            del_words.append(word.rstrip())
    f.close()

    # deleting titles with unwanted words
    for j in del_words:
        os.system(f'sed -i "/\<{j}\>/Id" ../DUFs/{dufID}/{pfamID}_all_titles_tmp.txt')

    # Retrieving hits from temporary copy of hit->titles file using words from word_dictionary created earlier as input
    try:
        with open(f"../DUFs/{dufID}/{pfamID}_putative_hits.txt", "w") as res:
            for j in range(len(hit_words)):
                for i in range(len(hit_words) - 1):
                    cmd1 = f'grep -i "{hit_words[j]}.*{hit_words[j - i - 1]}" ../DUFs/{dufID}/{pfamID}_all_titles_tmp.txt'
                    os.system(cmd1)
                    putative_hits = os.popen(cmd1).read()
                    res.write(putative_hits.strip() + '\n')
        res.close()

        # Deleting any empty lines from file containing all putative hits
        os.system(f"sed -i /^$/d ../DUFs/{dufID}/{pfamID}_putative_hits.txt")

        # deleting temporary copy of hit->titles file and file containing words deleted to generate filtered word dictionary
        os.remove(f"../DUFs/{dufID}/{pfamID}_all_titles_tmp.txt")
        os.remove(f"../DUFs/{dufID}/{pfamID}_del_list.txt")

        status = 'Putative annotations obtained'

    except:
        status = "Error during annotation"

    return status

# blastHitRst = sys.argv[1]
# annotate(blastHitRst)