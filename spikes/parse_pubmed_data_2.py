#!/usr/bin/env python
# encoding: utf-8
import re

from Bio import Entrez,Medline

#get all abstracts
txt=open("TGF-β - IL23r.txt","r")
records=Medline.parse(txt)
#title,firstAuthor,date,pubmedID
for r in records:
	print r["AB"]