import nltk
_POS_TAGGER = 'taggers/maxent_treebank_pos_tagger/english.pickle' 

posTagger=nltk.data.load(_POS_TAGGER)


abstract="""T cell functional differentiation is mediated by lineage-specific transcription factors. T helper 17 (Th17) has been recently identified as a distinct Th lineage mediating tissue inflammation. Retinoic acid receptor-related orphan receptor gamma (ROR gamma) was shown to regulate Th17 differentiation; ROR gamma deficiency, however, did not completely abolish Th17 cytokine expression. Here, we report Th17 cells highly expressed another related nuclear receptor, ROR alpha, induced by transforming growth factor-beta and interleukin-6 (IL-6), which is dependent on signal transducer and activator of transcription 3. Overexpression of ROR alpha promoted Th17 differentiation, possibly through the conserved noncoding sequence 2 in Il17-Il17f locus. ROR alpha deficiency resulted in reduced IL-17 expression in vitro and in vivo. Furthermore, ROR alpha and ROR gamma coexpression synergistically led to greater Th17 differentiation. Double deficiencies in ROR alpha and ROR gamma globally impaired Th17 generation and completely protected mice against experimental autoimmune encephalomyelitis. Therefore, Th17 differentiation is directed by two lineage-specific nuclear receptors, ROR alpha and ROR gamma."""

abstracts=[
"""TGF-beta acts to upregulate IL-23R expression, thereby conferring responsiveness to IL-23""",
"""However, unexpectedly, the up-regulation of IL-23R and induction of Th17 differentiation by IL-6 and IL-23 were almost completely inhibited by anti-TGF-beta""",
"""These results suggest that the induction of IL-23R and Th17 differentiation by IL-6 and IL-23 is mediated through endogenously produced TGF-beta."""
]

# sents=abstract.split(".")
# for s in sents:
# 	ss=nltk.word_tokenize(s)
# 	print nltk.pos_tag(ss)

interestingPos=["J","V","MOD"]
def accepting_pos(pos):
	if pos[1] in interestingPos:
		return True
	for p in interestingPos:
		if pos[1].startswith(p):
			return True
	return False
	
tokenizedSentenceByKw1=[
"T helper 17 (Th17)" , "has" , "been" , "recently" , "identified" , "as" , "a" , "distinct" , "Th lineage" , "mediating" , "tissue inflammation"]
pos=nltk.pos_tag(tokenizedSentenceByKw1)
#extract verbs
verbs= [x for x in pos if accepting_pos(x)]


tokenizedSentenceByKw2=[" Retinoic acid receptor-related orphan receptor gamma" , "ROR gamma" , "was" , "shown" , "to" , "regulate" , "Th17 differentiation" , "ROR gamma" , "deficiency" , "however" , "did" , "not" , "completely" , "abolish" , "Th17" , "cytokine" , "expression"
]



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