import nltk
_POS_TAGGER = 'taggers/maxent_treebank_pos_tagger/english.pickle' 

posTagger=nltk.data.load(_POS_TAGGER)


abstracts=[
"""TGF-beta acts to upregulate IL-23R expression, thereby conferring responsiveness to IL-23""",
"""However, unexpectedly, the up-regulation of IL-23R and induction of Th17 differentiation by IL-6 and IL-23 were almost completely inhibited by anti-TGF-beta""",
"""These results suggest that the induction of IL-23R and Th17 differentiation by IL-6 and IL-23 is mediated through endogenously produced TGF-beta."""
]


interestingPos=["J","V","MOD","CELL","SPECIE","TF","RECEPTOR"]
def accepting_pos(pos):
	if pos[1] in interestingPos:
		return True
	for p in interestingPos:
		if pos[1].startswith(p):
			return True
	return False
	

#tagger to recognize NT from the domain
ntr={
"T helper 17 (Th17)":"CELL",
"tissue inflammation":"COND",
'ROR gamma':"TF",
'Th17 differentiation':"PROCESS",
"IL-23":"SPECIE",
"TGF-beta":"SPECIE",
"IL-23R":"RECEPTOR",
' Retinoic acid receptor-related orphan receptor gamma':"TF",
}
baseline_tagger=nltk.UnigramTagger(model=ntr,backoff=posTagger)

for a in abstracts:
	#test on a sentence about tgfb and IL23R
	text = nltk.word_tokenize(a)
	
	tagged=baseline_tagger.tag(text)
	print tagged
	print [x for x in tagged if accepting_pos(x)]