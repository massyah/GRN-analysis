from gensim import corpora,models,similarities
import psycopg2
import logging
import numpy,datetime,random,re
import postgres_indexing as pg


word_b=re.compile(r'\w+',re.U)
stoplist=set(load_latest_data_named("stop_word_english_cornell"))

#some terms to recognize
terms=list(pg.load_latest_data_named("th17_human_terms"))
terms_only_chemicals=[x for x in terms if x.isChemical]
text_to_terms=build_term_dict(terms)
termsRE=build_re(text_to_terms.keys())



def load_lsi_data():
	#loaded in 7s
	global lsi,dictionary,corpus,index,pmids,dates,titles
	lsi,dictionary,corpus,index,pmids,dates,titles=cPickle.load(open("gensim_temp.dat"))
	lsi.projection.u=lsi.projection.u.copy('F')
	
def split_word_boundary(txt):
	txt=txt.replace("-","_") #should be corrected by NTR but let's keep it for testing
	words= word_b.findall(txt)
	return [w for w in words if w not in stoplist]
	
def build_lsi_for_abstracts(numAbstract=None):
	global dictionary,lsi,index,pmids,titles,documents,dates,corpus
	conn=psycopg2.connect("dbname=th17 user=hayssam password=mjnfqto")
	cur=conn.cursor()

	if numAbstract:
		q="""SELECT pmid,title,da,abstract from "Publication" LIMIT %d;"""%(numAbstract)
	else:
		q="""SELECT pmid,title,da,abstract from "Publication";"""
	cur.execute(q)
	documents=[]
	titles=[]
	pmids=[]
	dates=[]
	for r in cur.fetchall():
		titles.append(r[1])
		pmids.append(r[0])
		dates.append(r[2])
		documents.append(normalizeText(r[1]+"\n"+r[3]))
	print "Index pmids",len(pmids)
	print "Normalization finished, Will compute LSI"
	sys.stdout.flush()
	texts = [split_word_boundary(document.lower()) for document in documents]

	# # remove words that appear only once
	# allTokens = sum(texts, [])
	# tokensOnce = set(word for word in set(allTokens) if allTokens.count(word) == 1)
	# texts = [[word for word in text if word not in tokensOnce] for text in texts]

	#build the lsi model
	dictionary=corpora.Dictionary(texts)
	corpus=[dictionary.doc2bow(text) for text in texts]
	lsi = models.LsiModel(corpus, id2word=dictionary.id2word, numTopics=200)

	 # transform corpus to LSI space and index it
	index = similarities.MatrixSimilarity(lsi[corpus])

def search(query):
	global dictionary,lsi,index,pmids,dates,titles
	# query="il-1 TNF-alpha"
	query=normalizeText(query)
	vec_bow=dictionary.doc2bow(query.lower().split())

	vec_lsi=lsi[vec_bow]
	sim=index[vec_lsi]
	res=map(lambda x:(x[1],pmids[x[0]],dates[x[0]],titles[x[0]]),list(enumerate(sim)))
	res.sort(key=lambda item:item[0],reverse=True)
	#count number of items with more than 0.5 sim score
	print len([x for x in res if x[0]>0.5]),"relevant docs"
	print "\n".join(map(str,res))

def canonical_name(term):
	return term.name.title().replace(" ","")

def replace_with_canonical(match):
	t=match.group(0).lower()
	if t not in text_to_terms:
		t="_"+t
	return canonical_name(text_to_terms[t])

def normalizeText(text):
	return termsRE.sub(replace_with_canonical,text)


def ranks_of_articles_for_query(query,articles,maxDate=None):
	global dictionary,lsi,index
	# query="il-1 TNF-alpha"
	query=normalizeText(query)
	vec_bow=dictionary.doc2bow(query.lower().split())

	vec_lsi=lsi[vec_bow]
	sim=index[vec_lsi]
	res=map(lambda x:(x[1],pmids[x[0]],dates[x[0]],titles[x[0]]),list(enumerate(sim)))
	# print "# NDOCS:",len(res)
	if maxDate!=None:
		res=[x for x in res if x[2]<=maxDate]
		print "#filtered",len(res)
	res.sort(key=lambda item:item[0],reverse=True)
	sorted_pmids=[x[1] for x in res]
	sorted_score=[x[0] for x in res]
	pmid_to_rank={}
	pmid_to_score=dict([(x[1],x[0]) for x in res])
	for p in articles:
		if p not in sorted_pmids:
			pmid_to_rank[p]=-1
			pmid_to_score[p]=-1
		else:
			pmid_to_rank[p]=sorted_pmids.index(p)
	return pmid_to_rank,pmid_to_score

def test_mednoza_cohesion(q):
	#remove all articles published after mendoza article
	pmid_cited_by_mendoza=[12461566,12547678,8700208,11298330,12094404,10760800,12778462,11516335,7796298,12947222,12797537,9120387,10190904,11015439,11457889,9548487,9590262,10706670,10605012,12193719,11739558,11283251,10851054,10097134,11752460,12242343,15261525,11397944,10679398,11907070,10725704,9200456,10201892,11175814,9531298,10358772,12500979,7612223,10761931,12006974,10661403,11374299,12479817,12648458,11714753,10820262];
	if q=="":
		q=load_latest_data_named("mendoza_2006_model_paragraph").decode("utf-8")
	mendoza_da=datetime.date(2006,4,17)
	# mendoza_da=None
	res,scores=ranks_of_articles_for_query(q,pmid_cited_by_mendoza,mendoza_da)
	scores=dict([(p,scores[p]) for p in pmid_cited_by_mendoza])
	# print scores
	print "Rank AVG STD",numpy.average([x[1] for x in res.items()]),numpy.std([x[1] for x in res.items()])
	print "Score AVG STD",numpy.average([x[1] for x in scores.items()]),numpy.std([x[1] for x in scores.items()])
	return sorted(res.items(),key=itemgetter(1))

def consistency_of_pmids(pmidsSet):
	q="\n".join([get_abstract_from_db(x) for x in pmidsSet])
	res,scores=ranks_of_articles_for_query(q,pmidsSet)
	scores=dict([(p,scores[p]) for p in pmidsSet])
	# print scores
	print "Rank AVG STD",numpy.average([x[1] for x in res.items()]),numpy.std([x[1] for x in res.items()])
	print "Score AVG STD",numpy.average([x[1] for x in scores.items()]),numpy.std([x[1] for x in scores.items()])
	return sorted(res.items(),key=itemgetter(1))
	
def consistency_of_pmids_by_size(pmidsSet):
	q="\n".join([get_abstract_from_db(x) for x in pmidsSet])
	res,scores=ranks_of_articles_for_query(q,pmidsSet)
	scores=dict([(p,scores[p]) for p in pmidsSet])
	# print scores
	avg_rank=numpy.average([x[1] for x in res.items()])
	avg_std=numpy.std([x[1] for x in res.items()])
	return len(pmidsSet),avg_rank,avg_std


corpus=None
dictionary=None
lsi=None
index=None
# load_lsi_data()
# pmids=[]
# titles=[]
# build_lsi_for_abstracts()# documents=[]
