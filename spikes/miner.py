#!/usr/bin/env python
# encoding: utf-8
import re,copy


test="In contrast, STAT-1 cannot be activated by IL-4 (Moriggl et al, 1998)"

allPublications={}
allSentences=[]
allTerms={}
allTermsRe=None
allPredicates=[]
uid=0

def citeKey_to_publication(citeKey):
	if citeKey.lower() in allPublications:
		return allPublications[citeKey.lower()]
	else:
		return None
		
class Publication(object):
	"""docstring for Publication"""
	def __init__(self, title,firstAuthor,date):
		global allPublications,uid
		super(Publication, self).__init__()
		self.uid=uid
		uid+=1
		self.title = title
		self.firstAuthor = firstAuthor
		self.date = date
		self.citeKey=str(date)+" "+firstAuthor.lower()
		if self.citeKey in allPublications:
			print "Publication key conflict for",self.citeKey
		allPublications[self.citeKey]=self
	def __str__(self):
		return "%s (%d)"%(self.citeKey,self.uid)
		
		
def paragraph_to_sentences(pub,paragraph):
	for s in paragraph.split("\n"):
		if len(s)<=1:
			continue
		Sentence(pub,s)
	
	
def highlighter_call_back(matchobj):
	return "<B>%s<B>"%(matchobj.group(0))
	
class Sentence(object):
	"""docstring for Sentence"""
	def __init__(self, publication,string):
		global allSentences,uid
		super(Sentence, self).__init__()
		self.string = string
		self.publication=publication
		self.terms=[]
		self.uid=uid
		uid+=1
		allSentences.append(self)
	def term_recognition(self):
		for m in allTermsRe.findall(self.string):
			self.terms.append(allTerms[m.lower()])
	def highlight(self):
		outputString=copy.copy(self.string)
		return str(self.uid)+":\t"+allTermsRe.sub(highlighter_call_back,outputString)
		


def sentences_with_terms(terms):
	results=[]
	setOtTerms=set(terms)	
	for sent in allSentences:
		setOfTermsInSent=set(sent.terms)
		if setOtTerms<=setOfTermsInSent:
			results.append(sent)
	return results
	
	
def _updateRE():
	global allTerms,allTermsRe
	keys = sorted(allTerms.keys(), key=len)
	keys.reverse()
	expression = []
	for item in keys:
		if item[0]=="_": #then word boundary are important
			expression.append("\\b"+re.escape(item[1:])+"\\b")
		else:
			expression.append(re.escape(item))
	allTermsRe = re.compile("(%s)" % "|".join(expression),re.IGNORECASE)
	
class Term(object):
	"""docstring for Term"""
	def __init__(self, name,alternativeNames):
		global allTerms,uid,allTermsRe
		super(Term, self).__init__()
		self.name = name
		self.alternativeNames = alternativeNames
		self.uid=uid
		uid+=1
		if name.lower() in allTerms:
			print "Name conflict for term",name
		allTerms[name.lower()]=self
		for an in alternativeNames:
			if an.lower() in allTerms:
				print "Name conflict for term",name,"already present in",allTerms
				assert(False)
			allTerms[an.lower()]=self
		_updateRE()
	
	def __eq__(self,other):
		return self.name==other.name

class Predicate(object):
	"""docstring for Predicate"""
	def __init__(self, name,obj,subj):
		global allTerms, uid
		super(Predicate, self).__init__()
		self.name = name
		self.obj = obj
		self.subj = subj
		self.uid=uid
		uid+=1
		allTerms.append(self)



mine=Publication("My findings","hayssam",2010)
moriggl=Publication("Activation of STAT proteins and cytokine genes in human Th1 and Th2 cells generated in the absence of IL-12 and IL-4","Moriggl",1998)

#interleukines
il4=Term("IL-4",["IL4"])
il6=Term("IL-6",["IL6"])

stat1=Term("STAT-1",["STAT1"])
stat1b=Term(u"STAT-1ß",["STAT1 beta",u"STAT1ß"])
stat3=Term("STAT-3",["STAT3"])
stat6=Term("STAT-6",["STAT6"])
th2=Term("Th2",["Th 2"])
ifng=Term(u"IFN-γ",["Interferon gamma","ifng"])
ifnga=Term(u"IFN-γa",["Interferon gamma a"])
il12=Term("IL-12",["IL12"])
rorgt=Term(u"ROR-γt",[u"RORγt","ROR gt"])
rora=Term(u"RORα",[u"RORa"])
wt=Term("wild type",["_WT"])
ahrko=Term("ahr-deficient",[])
il2ko=Term("il2-/-",[])
socs1=Term("SOCS-1",["SOCS1","SSI-1","STAT-induced STAT inhibitor 1","Suppressor of cytokine signaling 1"])
stat1ko=Term("STAT1-deficient",["STAT1-/-"])
ige=Term("IgE",[])

#cell types
cd4p=Term("CD4+",[])
th1=Term("Th1",[])
#Receptors
#TODO:Make a template for STATs, and receptors
il23r=Term("IL23r",["IL-23r","IL23-r","IL-23 r","IL-23-r"])
ifngr=Term(u"IFNγR",[u"IFNγ receptor"])


Sentence(moriggl,u"In contrast, STAT-1 cannot be activated by IL-4.")
Sentence(moriggl,u"""
However, STAT1 and STAT3 can be activated by IL-12 and IFN-γ (21), two major players in the induction of Th1 immunity (18, 20).""")
Sentence(moriggl,"""STAT1 and STAT3 are not activated by IL-4 (21), the cytokine driving Th2 cell differentiation.""")

pub_venkat=Publication("Repression of IL-4-induced gene expression by IFN-gamma requires Stat1 activation","Venkataraman",1999)


paragraph_to_sentences(pub_venkat,u"""
IFN-γ antagonizes many physiological responses mediated by IL-4, including the inhibition of IL-4-induced IgE production.
This event is largely mediated at the level of transcription.
We observed that the IL-4 response element of the germline epsilon promoter is sufficient to confer IFN-γ-mediated repression onto a reporter construct.
The inhibitory effects were observed in both lymphoid and nonlymphoid cell lines.
Stat1, which is activated by IFN-γ, cannot recognize the Stat6-specific IL-4 response element in the promoter.
Hence, competitive DNA binding does not seem to be the underlying mechanism for the inhibitory effect.
This is supported by the observation that inhibition is not seen at early time points, but requires prolonged IFN-γ treatment.
IFN-γ stimulation results in a loss of IL-4-induced Stat6 tyrosine phosphorylation, nuclear translocation, and DNA binding.
Using the fibrosarcoma cell line U3A, which lacks Stat1, we demonstrated that the transcription activation function of Stat1 is required for the IFN-γ-mediated repression.
Repression was restored by overexpression of Stat1, but not Stat1ß, in U3A cells.
Treatment with IFN-γ, but not IL-4, specifically up-regulates the expression of SOCS-1 (silencer of cytokine signaling), a recently characterized inhibitor of cytokine signaling pathways, such as IL-6 and IFN-γ.
Overexpression of SOCS-1 effectively blocks IL-4-induced Stat6 phosphorylation and transcription.
This suggests that IFN-γ-mediated repression of IL-4-induced transcription is at least in part mediated by SOCS-1.""")


pub_saito=Publication(u"IFN regulatory factor-1-mediated transcriptional activation of mouse STAT-induced STAT inhibitor-1 gene promoter by IFN-γ","Saito",2000)
paragraph_to_sentences(pub_saito,u"""
We also showed that IFN-γ-induced SSI-1 mRNA more strongly than IL-6 in NIH-3T3 fibroblasts and that this IFN-γ effect was mediated by Stat1.
As previously shown that SSI-1 mRNA was induced by IL-6 in a Stat3-dependent manner in M1 cells and that activation of Stat1 was a main pathway of IFN-γ signaling, we attempted to verify whether Stat1 mediates this SSI-1 mRNA induction by IFN-γ.
Responsiveness to IFN-γ of embryonal fibroblasts prepared from Stat1-deficient mice was analyzed.

As shown in Fig. 1D, Stat1-/- fibroblasts did not respond to IFN-γ stimulation in terms of SSI-1 mRNA induction, and concordantly, IRF-1 mRNA induction was diminished, while Stat1+/+ fibroblasts responded to IFN-γ as expected.
These results indicated that IFN-γ induced SSI-1 via a Stat1-dependent pathway.
""")



pub_diehl=Publication(u"Inhibition of Th1 differentiation by IL-6 is mediated by SOCS1","Diehl",2000)
paragraph_to_sentences(pub_diehl,u"""
While IL-6-directed CD4+ Th2 differentiation is mediated by IL-4, inhibition of Th1 differentiation by IL-6 is independent of IL-4.
IL-6 upregulates suppressor of cytokine signaling 1 (SOCS1) expression in activated CD4+ T cells, thereby interfering with signal transducer and activator of transcription 1 (STAT1) phosphorylation induced by interferon γ (IFNγ).
Inhibition of IFNγ receptor-mediated signals by IL-6 prevents autoregulation of IFNγ gene expression by IFNγ during CD4+ T cell activation, thereby preventing Th1 differentiation.
Thus, IL-6 promotes CD4+ Th2 differentiation and inhibits Th1 differentiation by two independent molecular mechanisms.
Although the role of IFNγ in Th1 differentiation is somewhat controversial, several studies have shown the requirement for IFNγ for complete Th1 differentiation ([8, 51, 50 and 55]).
Thus, it was possible that the low levels of IFNγ that were produced during CD4+ T cell differentiation in the presence of IL-6 were responsible for impaired Th1 differentiation.
To address this possibility, CD4+ T cells were stimulated in the presence or absence of IL-6 and exogenous IFNγ, washed, and restimulated.
Low levels of IFNγ were produced by effector cells differentiated with IL-6, even if exogenous IFNγ was present during differentiation (Figure 4A).
Thus, the presence of exogenous IFNγ during stimulation could not compensate for the inhibitory effect on Th1 differentiation exerted by IL-6.
However, these results would be expected if the negative regulatory effect of IL-6 on Th1 differentiation was due to an impairment of the IFNγ receptor (IFNγR)-signaling pathways involved in the upregulation of IFNγ gene expression.
IL-6 Inhibits Signaling through the IFNγ Receptor by Upregulating SOCS1
IL-6 Fails to Inhibit IFNγ Production in SOCS1-Deficient Mice
IL-6 Interferes with IFNγ Receptor-Mediated IFNγ Production
IL-6 Upregulates SOCS1 Expression and Inhibits IFNγ-Induced STAT1 Phosphorylation
To show that elevated levels of SOCS1 induced by IL-6 in CD4+ T cells could inhibit IFNγR-mediated signals, we examined Tyr- phosphorylation of STAT1.
CD4+ T cells were stimulated for 24 hr in the presence or absence of IL-6.
No phosphorylated STAT1 was detected in unstimulated CD4+ T cells, but significant levels were found in activated cells (Figure 6E).
This phosphorylation was primarily induced by endogenous IFNγ since the presence an anti-IFNγ mAb blocked STAT1 phosphorylation ( Figure 6E).
Cells activated in the presence of IL-6 also exhibited decreased levels of phospho-STAT1 ( Figure 6E).
To better examine the inhibition of IFNγ-induced STAT1 phosphorylation by IL-6 activated CD4+, T cells were rested to decrease phospho-STAT1 background and then treated with IFNγ.
STAT1 phosphorylation was rapidly induced by IFNγ in cells stimulated in the absence of IL-6 while no significant induction of STAT1 phosphorylation was induced in cells activated in the presence of IL-6 (Figure 6F).
The levels of STAT1 protein were similar in all the conditions.
Thus, the IFNγ-signaling pathway is impaired in CD4+ T cells differentiated in the presence of IL-6.
To further demonstrate that IL-6 inhibited IFNγ production by CD4+ T cells by enhancing SOCS1 expression, we examined whether IL-6 failed to inhibit IFNγ production in SOCS1-deficient Th cells (SOCS1−/−) ([23]).
CD4+ T cells from wild-type or SOCS1−/− littermates were stimulated in the presence or absence of IL-6.
The presence of IL-6 completely blocked the early production of IFNγ by wild-type CD4+ T cells, but IL-6 did not affect IFNγ production by CD4+ T cells that lack SOCS1 (Figure 7A).
Similar results were obtained when cells were differentiated in the presence of IL-12 ( Figure 7B).
These data demonstrate that IL-6 inhibited IFNγ signaling and IFNγ gene expression by upregulating SOCS1 gene expression.
""")

pub_kotenko=Publication(u"Jak-Stat signal transduction pathway through the eyes of cytokine class II receptor complexes","Kotenko",2000)
paragraph_to_sentences(pub_kotenko,
"IFN-γ triggers Stat1 and, in some cases, dependent on the cell type, specific conditions and/or ratio of receptor subunits expressed, a small amount of Stat3.")

pub_goodbourn=Publication("","Goodbourn",2000)
paragraph_to_sentences(pub_goodbourn,
u"IFN-a/β activates the STAT1/STAT2 Phosphorylated dimers, which complexes later with the DNA binding protein p58 to form an heterotrimeric complex called ISGF3, which binds to the promoter region of IFN-a/β genes.")

for s in allSentences:
	s.term_recognition()

for s in sentences_with_terms([ifngr]):
	print s.highlight()


#verify that wt in a subword is not considered as a term
#display all sentences with term highlighted
