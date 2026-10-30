import urllib.request
url = "https://rest.uniprot.org/uniprotkb/stream?query=organism_id:9606%20AND%20reviewed:true&format=fasta"
data = urllib.request.urlopen(url, timeout=300).read()
open("data/processed/_uniprot_human_reviewed.fasta", "wb").write(data)
print("saved", len(data))
