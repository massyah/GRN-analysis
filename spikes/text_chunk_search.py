#!/usr/bin/env python
# encoding: utf-8

import re,copy
from appscript import *
txt=open("python_re.txt","r").read().decode("utf-8")
dt=app("DEVONThink Pro")

def shortest_chunk():
	query=re.compile("python(.*?)regular expression(.*)special(.*)whitespace",re.IGNORECASE|re.DOTALL)
	shortestSpan=float('inf')
	match=None
	for m in query.findall(txt):
		if len(m)<=shortestSpan:
			shortestSpan=len(m[0])+len(m[1])
			match=m
	print shortestSpan
	# print match

#n-grams approach
kw=["python","regular expression","special","matches","whitespace","character"]
kw_to_index={}
i=0
for k in kw:
	kw_to_index[k]=i
	i+=1
	
def tokenize():
	kwtxt=[]
	for k in kw:
		r=re.compile(k,re.IGNORECASE)
		for m in r.finditer(txt):
			kwtxt.append((k,m.span()[0]))
	kwtxt.sort(key=lambda x:x[1])
	kwtxt2=[]
	i=0
	for k in kwtxt:
		kwtxt2.append((i,k[0],k[1]))
		i+=1
	
def view_sentence(pos):
	print txt[pos-300:pos+300]
	
def view_context(chunk_id):
	view_sentence(kwtxt2[chunk_id][2])
	
def distance_between_kw(k1,k2):
	chunks=[x for x in kwtxt2 if x[1] in [k1,k2]]
	minDist=float('inf')
	for i in range(len(chunks)-1):
		if chunks[i][1]==chunks[i+1][1]:
			continue
		dist=chunks[i+1][2]-chunks[i][2]
		if dist<minDist:
			minDist=dist
			span=(chunks[i][2],chunks[i+1][2])
	# print txt[span[0]-200:span[1]+200]
	print k1,k2,span,minDist

#to pick the two keywords that are close to each other
def all_distances():
	distances=[]
	for i in range(len(kw)):
		distances.append([('inf',0)]*len(kw))

	for i in range(len(kwtxt2)-1):
		k1=kwtxt2[i][1]
		k2=kwtxt2[i+1][1]
		row=min(kw_to_index[k1],kw_to_index[k2])
		col=max(kw_to_index[k1],kw_to_index[k2])
		dist=kwtxt2[i+1][2]-kwtxt2[i][2]
		if distances[row][col][0]>dist:
			distances[row][col]=(dist,kwtxt2[i][2])
	print ("\t"*2).join([str(x)[:6] for x in kw])
	for r in distances:
		print ("\t"*2).join([str(x[0]) for x in r])
	return distances
#grouping the chunks by a 600 char window
def chunk_groups():
	groupedChunks=[]
	for i in range(len(kwtxt2)):
		thisChunk=[kwtxt2[i]]
		anchorPos=kwtxt2[i][2]
		imin=i-1
		iplus=i+1
		while imin>0 and anchorPos-kwtxt2[imin][2]<=300:
			thisChunk.insert(0,kwtxt2[imin])
			imin-=1
		while iplus<len(kwtxt2) and kwtxt2[iplus][2]-anchorPos<=300:
			thisChunk.append(kwtxt2[iplus])
			iplus+=1
		if len(thisChunk)>1:
			groupedChunks.append((anchorPos,thisChunk))
	return groupedChunks

def paragraph_scoring(txt):
	pars=txt.split("\n")
	for i in range(1,len(pars)-1):
		kwtxt=[]
		withNeigh=pars[i-1]+pars[i]+pars[i+1]
		#search for kw in neigh
		for k in kw:
			r=re.compile(k,re.IGNORECASE)
			for m in r.finditer(withNeigh):
				kwtxt.append((k,m.span()[0]))
		kwtxt.sort(key=lambda x:x[1])
		if len(kwtxt)>0:
			print kwtxt
			print withNeigh[:90],"....",withNeigh[-90:]
def topic_sentences(txt):
	pars=txt.split("\n")
	tsents=[]
	for p in pars:
		tsents.append(p.split(".")[0])
	tsents.sort()
	return tsents
	
contextSize=90
def devon_search():
	tsents=[]
	query=["python","re","whitespace"]
	for db in dt.databases():
		rec=dt.search(" ".join(query),in_=db.root())
		print len(rec)
		for r in rec:
			if r.word_count()<=20000:
				print r.name()
				txt=r.plain_text()
				for k in query:
					reg=re.compile(k,re.IGNORECASE)
					print k,1.0*len(reg.findall(txt))/r.word_count()
				# tsents.extend(topic_sentences(txt))
	# print tsents
	# print "\n".join([x for x in tsents if  "Basic" in x])

#the last keyword must appear
#the previous keyword are just here to compute anchors

#order the keywords by the number of matches
#start to present paragraphs with the kw with lowest occurence
#assuming that some keywords are here to find the document, while the others are there to find a specific sentence, paragraph

#try to find first in the topic sentences then in the other sentences

# def score(g):
# 	variety=len(set([x[1] for x in g[1]]))
# 	return 1.0*sum([kw_to_index[x[1]] for x in g[1]])/len(g[1])*variety
# groups=chunk_groups()
# groups.sort(key=lambda x:score(x),reverse=True)
# for g in groups[:30]:
# 	print score(g),g

# paragraph_scoring()