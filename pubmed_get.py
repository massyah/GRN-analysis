#!/usr/bin/env python
# encoding: utf-8

from Bio import Entrez,Medline
def get_abstract(query,file_name,previewOnly=False):
	Entrez.email = "massyah@gmail.com"
	search_results = Entrez.read(Entrez.esearch(db="pubmed",
																							term=query,
																							reldate=2000, datetype="pdat",
																							usehistory="y"))
	count = int(search_results["Count"])
	print "Found %i results" % count
	if previewOnly:
		return
	batch_size = 10
	out_handle = open(file_name+".txt", "w")
	for start in range(0,count,batch_size):
			end = min(count, start+batch_size)
			print "Going to download record %i to %i" % (start+1, end)
			fetch_handle = Entrez.efetch(db="pubmed",
																	 rettype="medline", retmode="text",
																	 retstart=start, retmax=batch_size,
																	 webenv=search_results["WebEnv"],
																	 query_key=search_results["QueryKey"])
			data = fetch_handle.read()
			fetch_handle.close()
			out_handle.write(data)
	out_handle.close()


def build_query(terms):
	termsQuery=[]
	for t in terms:
		termList=t.alternativeNames+[t.name]
		scope=["%s [Title/Abstract]"%(term) for term in termList]
		termsQuery.append("("+" OR ".join(scope)+")")
	return " AND ".join(termsQuery)
	
def build_file_name(terms):
	return "entrez "+" - ".join([t.name for t in terms])
	
def pubmed_get(terms):
	get_abstract(build_query(terms), build_file_name(terms))
def pubmed_preview(terms):
	q=build_query(terms)
	print q
	get_abstract(build_query(terms),None,previewOnly=True)