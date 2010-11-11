#!/usr/bin/env python
# encoding: utf-8


from whoosh.index import create_in,open_dir
from whoosh.fields import *
from whoosh.searching import *
from whoosh.qparser import QueryParser,SimpleParser

from Bio import Entrez,Medline
from whoosh.qparser import MultifieldParser

import os

project="Th17 differentiation"
termFile=project+"/articles_th17.txt"
pubmed_files=[f for f in os.listdir("./"+project) if f.startswith("entrez ")]
# pubmed_files=["entrez 16648837.txt"]
# print [type(x.decode("utf-8")) for x in pubmed_files]

schema = Schema(
				title=TEXT(stored=True), 
				path=ID(stored=True), 
				abstract=TEXT(stored=True),
				authors=TEXT(stored=True),
				pmid=TEXT(stored=True),
				dateAdded=DATETIME(stored=True)
			)



# writer = ix.writer()
# writer.add_document(title=u"First document", path=u"/a", content=u"This is the first document we've added!")
# writer.add_document(title=u"Second document", path=u"/b",content=u"The second one is even more interesting!")
# writer.commit()

def clean_index():
	ix = create_in("indexdir", schema)
def index():
	ix = open_dir("indexdir")
	writer = ix.writer()
	for pfile in pubmed_files:
		print "parsing",pfile
		txt=open(project+"/"+pfile,"r")
		records=Medline.parse(txt)
		for r in records:
			if "AB" not in r:
				continue
			authors=""
			if "FAU" in r:
				authors+=",".join(r["FAU"])
			elif "AU" in r:
				authors+=",".join(r["AU"])
			else:
				firstAuthor="Unknown"
			date=datetime.datetime.strptime(r["DA"],"%Y%m%d")
			title=r["TI"]
			pmid=r["PMID"].decode("utf-8")

			writer.add_document(
				title=title.decode("utf-8"),
				path=pfile.decode("utf-8"),
				abstract=r['AB'].decode("utf-8"),
				authors=authors.decode("utf-8"),
				pmid=pmid,
				dateAdded=date
				)
	writer.commit()
	print "Index contain",ix.doc_count()


def search(q):
	print "Index contain",ix.doc_count()
	searcher = ix.searcher()
	mparser = MultifieldParser(["title", "abstract","authors","pmid"], schema=schema)
	# query = QueryParser("abstract").parse("aberrant cytokine mangan")
	# query = QueryParser("authors").parse("mangan")
	query=mparser.parse(q)
	# query=QueryParser("abstract").parse(q)

	# query=mparser.parse(u"\"IL-1β\"")

	results = searcher.search(query)

	resItems=results.items()
	print len(resItems),"doc found"
	for i in range(len(resItems)):
		print results[i]["pmid"],results[i]["dateAdded"],results[i]["title"],resItems[i][1]
	


def add_sample():
	ix = open_dir("indexdir")
	writer = ix.writer()
	abs1=u"Everything is XXXXX. This is about IL-12."
	abs2=u"Everything is not XXXXX so anymore near. but XXXXX can it be? and XXXXX again? This is about IL-21 and IL-12 and so much more."

	writer.add_document(abstract=abs1,title=u"abs1",pmid=1)
	writer.add_document(abstract=abs2,title=u"abs2",pmid=2)
	writer.commit()
	print "Index contain",ix.doc_count()

def pmid_to_keywords(pmidval):
	searcher=ix.searcher()
	docn=searcher.document_number(pmid=pmidval)
	print list(searcher.key_terms([docn], "abstract"))

# pmid_to_keywords(16648837)
# search(u"\"TGF-β\"")


#test unicode OK
#test ranking OK
#test if sentence proximity affect ranking NO
#make compatible with synonym expansion
#think about a sentence based indexing
search(u"\"TGF-β\"")