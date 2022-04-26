import urllib.request

f=open("Apoptosis Network (SBML).sbml")

lines=f.readlines()

f.close()
links=[]
for line in lines:
    if "qual:name" in line:
        fst = line.index("qual:name") + 10
        temp = line[fst:]
        name = line[fst+1:-3]
        print(name)
    if "UniProt ID:" in line:
        fst = line.index("UniProt ID:") + 10
        temp = line[fst:]
        lst = temp.index('"')
        ID = line[fst:lst + fst]
        print(ID)
    if "href" in line:
        fst=line.index("href")+6
        temp=line[fst:]
        lst=temp.index('"')
        link=line[fst:lst+fst]+".fasta"

        if "uniprot" in link:
            r=urllib.request.urlopen(link)
            links.append(link)
            content=r.read()
            content=str(content)
            content=content.split("\\n")
            #print(content)
            sequence="".join(content[1:-1])
            print(sequence)

print(len(links))