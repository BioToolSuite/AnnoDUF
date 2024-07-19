import os
import sys

def annotate(dufID, pfamID, c):
    putativeHitsFile = f"../DUFs/{dufID}/{pfamID}_Split_{c}_putative_hits.txt"
    blastHitRstFile = f"../DUFs/{dufID}/{pfamID}_Split_{c}_all_titles.txt"
    allTempTitlesFile = f"../DUFs/{dufID}/{pfamID}_Split_all_titles_tmp.txt"
    wordDictFile = f"../DUFs/{dufID}/{pfamID}_Split_{c}_word_dict.txt"
    delListFile = f"../DUFs/{dufID}/{pfamID}_Split_{c}_del_list.txt"
    if c == 0:
        putativeHitsFile = f"../DUFs/{dufID}/{pfamID}_putative_hits.txt"
        blastHitRstFile = f"../DUFs/{dufID}/{pfamID}_all_titles.txt"
        allTempTitlesFile = f"../DUFs/{dufID}/{pfamID}_all_titles_tmp.txt"
        wordDictFile = f"../DUFs/{dufID}/{pfamID}_word_dict.txt"
        delListFile = f"../DUFs/{dufID}/{pfamID}_del_list.txt"
    
    # print("**Annotation under process**")
    # creating a temporary copy of hit->titles file
    cmd = f'sort -f {blastHitRstFile} | uniq -i > {allTempTitlesFile}'
    os.system(cmd)

    hit_words = []
    del_words = []

    # appending words from filtered word dictionary
    with open(wordDictFile) as f:
        for word in f:
            hit_words.append((word.rstrip().split("'"))[1])
    f.close()

    # appending words which were deleted to generate filtered word dictionary
    with open(delListFile) as f:
        for word in f:
            del_words.append(word.rstrip())
    f.close()

    # deleting titles with unwanted words
    for j in del_words:
        os.system(f'sed -i "/\<{j}\>/Id" {allTempTitlesFile}')

    # Retrieving hits from temporary copy of hit->titles file using words from word_dictionary created earlier as input
    try:
        with open(putativeHitsFile, "w") as res:
            for j in range(len(hit_words)):
                for i in range(len(hit_words) - 1):
                    cmd1 = f'grep -i "{hit_words[j]}.*{hit_words[j - i - 1]}" {allTempTitlesFile}'
                    os.system(cmd1)
                    putative_hits = os.popen(cmd1).read()
                    res.write(putative_hits.strip() + '\n')
        res.close()

        # Deleting any empty lines from file containing all putative hits
        os.system(f"sed -i /^$/d {putativeHitsFile}")

        # deleting temporary copy of hit->titles file and file containing words deleted to generate filtered word dictionary
        os.remove(allTempTitlesFile)
        os.remove(delListFile)

        status = 'Putative annotations obtained'

    except:
        status = "Error during annotation"

    return status
