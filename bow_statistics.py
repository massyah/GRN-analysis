bow_statistics={}
for sent in allSentences:
	terms=list(set(sent.terms))
	terms.sort(key=lambda x:x.name)
	terms=tuple(terms)
	if terms not in bow_statistics:
		bow_statistics[terms]=0
	bow_statistics[terms]+=1
	
bow_statistics=sorted(bow_statistics.iteritems(), key=lambda x:x[1], reverse=True)
for k,v in bow_statistics:
	print map(lambda x:x.name,k),"\t",v
	