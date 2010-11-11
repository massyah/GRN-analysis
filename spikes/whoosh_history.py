from appscript import *
from whoosh.index import create_in,open_dir
from whoosh.fields import *
from whoosh.searching import *
from whoosh.qparser import QueryParser,SimpleParser

from Bio import Entrez,Medline
from whoosh.qparser import MultifieldParser

import os



sa=app("Safari")
def safari_open_docs():
	docs=[]
	# for w in sa.windows():
	for w in [sa.windows()[0]]:
		if w.name()=="Downloads":
			continue
		for t in w.tabs():
			if "google" in t.URL():
				continue
			docs.append( {"title":t.name(), "URL":t.URL(), "content":t.text(), "dateAdded":datetime.datetime.now()
})
	return docs


schema = Schema(
				title=TEXT(stored=True), 
				content=TEXT(stored=True),
				url=TEXT(stored=True),
				dateAdded=DATETIME(stored=True)
			)

def clean_index():
	ix = create_in("history", schema)

def index_safari():
	ix = open_dir("history")
	writer = ix.writer()
	# try:
	for doc in safari_open_docs():
		print "parsing",doc["URL"]
		writer.add_document(
			title=doc["title"],
			url=doc["URL"].decode("utf-8"),
			dateAdded=doc["dateAdded"],
			content=doc["content"]
			)
	writer.commit()
	# except Exception as inst:
	# 	print "exception:",inst
	# 	writer.cancel()
	print "Index contain",ix.doc_count()
	
# clean_index()
# index_safari()

def search():
	while True:
		q=raw_input()
		print "search for ",q
		