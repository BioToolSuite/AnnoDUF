import requests
proxies = {'https':'http://proxy.ibab.ac.in:3128', 'http':'http://proxy.ibab.ac.in:3128'}
def uniprot(searchText):
    searchText = searchText.replace(" ", "%20")
    url = f'https://rest.uniprot.org/uniprotkb/search?query={searchText}'
    data = requests.get(url, proxies=proxies).json()
    if "messages" in data or len(data['results']) == 0:
        return "No hits!"
    else:
        return f"https://www.uniprot.org/uniprotkb?query={searchText}"
uniprot("DUF3255 family protein")