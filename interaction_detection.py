verbs="activat binds inhibit repress regulate suppress express associate block mediated inactivate contain induce modify overexpress promote stimulate enhanc".split(" ")
actors=[il_6,stat3]
actors_keys=[]
for a in actors:
	actors_keys.append(a.name)
	actors_keys.extend(a.alternativeNames)

# verbs.extend(target)
actor_re=build_re(actors_keys)
inter_re=build_re(verbs)

actorsButNoVerbs=[]
def score_text(t):
	this_actors=set()
	for m in actor_re.findall(t):
		this_actors.add(allTerms[m.lower()])

	if len(this_actors)>=len(actors):
		verbs=inter_re.findall(t)
		if len(verbs)>0:
			return 1,verbs
		else:
			actorsButNoVerbs.append(t)
	return 0,[]
verbs_to_sentences={}
#get sentences from the abstracts
conn=psycopg2.connect("dbname=th17 user=hayssam password=mjnfqto")
cur=conn.cursor()
#article with same PMID already here?
cur.execute("SELECT pmid,da,abstract FROM \"Publication\"" )
for r in cur.fetchall():
	ab=r[2]
	for sent in ab.split("."):
		score,verbs=score_text(sent)
		if score>0:
			for v in verbs:
				if v not in verbs_to_sentences:
					verbs_to_sentences[v]=[]
				verbs_to_sentences[v].append((r[1],r[0],sent))
		
print "-----------------Not recognized"
print "\n\t".join(actorsButNoVerbs)

for k,v in verbs_to_sentences.items():
	print "-"*8,"Verb:",k
	print "\n".join([x[2] for x in v])