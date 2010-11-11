from Bio import Entrez,Medline

txt=open("pmid 19379825.txt","r")
records=Medline.parse(txt)
for r in records:
	print r['FT']
