#!/usr/bin/env python
# encoding: utf-8
import re,copy
import omg_interface as omg

from IPython import ColorANSI
from IPython.genutils import Term
tc = ColorANSI.TermColors()



allPublications={}
allSentences=[]
allTerms={}
allTermsRe=None
allPredicates=[]
uid=0
sentencesUid=0


def citeKey_to_publication(citeKey):
	if citeKey.lower() in allPublications:
		return allPublications[citeKey.lower()]
	else:
		return None
		
class Publication(object):
	"""docstring for Publication"""
	def __init__(self, title,firstAuthor,date,pubmedID=None):
		global allPublications,uid
		super(Publication, self).__init__()
		self.uid=uid
		uid+=1
		self.title = title
		self.firstAuthor = firstAuthor
		self.date = date
		self.citeKey=str(date)+" "+firstAuthor.lower()
		self.pubmedID=pubmedID
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
def pyHighlighter_call_back(matchobj):
	return tc.Red+matchobj.group(0)+tc.Normal
	
class Sentence(object):
	"""docstring for Sentence"""
	def __init__(self, publication,string):
		global allSentences,uid,sentencesUid
		super(Sentence, self).__init__()
		self.string = string
		self.publication=publication
		self.terms=[]
		self.uid=sentencesUid
		sentencesUid+=1
		allSentences.append(self)
	def term_recognition(self):
		for m in allTermsRe.findall(self.string):
			self.terms.append(allTerms[m.lower()])
	def highlight(self,iPython=True):
		outputString=copy.copy(self.string)
		if iPython:
			return str(self.uid)+":\t"+allTermsRe.sub(pyHighlighter_call_back,outputString)
		else:
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
		self.omg=None
		self.alternativeNames = alternativeNames
		self.uid=uid
		uid+=1
		if name.lower() in allTerms:
			print "Name conflict for term",name
		allTerms[name.lower()]=self
		for an in alternativeNames:
			if an.lower() in allTerms:
				print "Name conflict for term",an,"already present in",allTerms.keys()
				assert(False)
			allTerms[an.lower()]=self
		_updateRE()
	
	def __eq__(self,other):
		return self.name==other.name

mine=Publication("My findings","hayssam",2010)

def parse_file():
	lines=open("articles_th17.txt").readlines()
	print len(lines)
	lastPub=None
	global allPublications,allTerms,allSentences,allPredicates,uid
	uid=0
	allPublications={}
	allTerms={}
	allSentences=[]
	allPredicates=[]
	line_number=1
	for l in lines:
		line_number+=1
		l=l.strip()
		if len(l)<1:
			continue
		l=unicode(l,"utf-8")
		if l.startswith("#t"):
			l=l[2:].strip()
			tt=l.split(",")
			alternatives=[]
			if len(tt)>1 and len(tt[1])>1:
				alternatives=tt[1:]
			tt[0]=tt[0].strip()
			for i in range(len(alternatives)):
				alternatives[i]=alternatives[i].strip()
				alternatives[i]=alternatives[i].strip("_")
			t=Term(tt[0],alternatives)
		elif l.startswith("#p"): #create a new pub
			l=l[2:].strip()
			lastPub=Publication(*l.split(","))
			print "Created",lastPub
		elif l.startswith("/"):
			continue
		else:
			assert lastPub,"No pub declared for line %d"%(line_number)
			if lastPub.pubmedID!=None: #then it is a pubmed abstract, we split it by sentences
				for li in l.split(". "):
					Sentence(lastPub,li)
			else:
				Sentence(lastPub,l)

			
	for s in allSentences:
		s.term_recognition()
	
def generate_variable_command(str):
	cmd=""
	for l in str.split("\n"):
		if len(l)<2:
			continue
		for v in l.split(","):
			v=v.strip()
			if len(v)<1:
				continue
			vn=v.replace("-","_")
			cmd+="%s=allTerms[\"%s\"]\n" %(vn,v)
	return cmd

def pprint_sentences(sents):
	for s in sents:
		print s.highlight()
	
def search(terms):
	sents=sentences_with_terms(terms)
	pprint_sentences(sents)
def publication_to_sentences(citeKey):
	pub=allPublications[citeKey]
	sents=[]
	for s in allSentences:
		if s.publication==pub:
			sents.append(s)
	pprint_sentences(sents)
	
parse_file()



def print_graph():
	for p in allPredicates:
		print p.__str__()

class Predicate(object):
	"""docstring for Predicate"""
	def __init__(self, name,obj,subj,direct=False,evidenceSentence=None):
		global allTerms, uid
		super(Predicate, self).__init__()
		self.name = name
		self.obj = obj
		self.subj = subj
		self.uid=uid
		self.direct=direct
		uid+=1
		allPredicates.append(self)
	def __str__(self):
		return "%s -> %s"%(self.obj.name,self.subj.name)





#first mappings
class Activation(Predicate):
	def __init__(self, obj,subj,direct=False,evidenceSentence=None):
		super(Activation, self).__init__("Activation",obj,subj,direct,evidenceSentence)
	def __str__(self):
		return "%s -+> %s"%(self.obj.name,self.subj.name)

class Inhibition(Predicate):
	def __init__(self, obj,subj,direct=False,evidenceSentence=None):
		super(Inhibition, self).__init__("Inhibition",obj,subj,direct,evidenceSentence)
	def __str__(self):
		return "%s -| %s"%(self.obj.name,self.subj.name)
	
def SignalingAxis(axis):
	for i in range(1,len(axis)):
		Activation(axis[i-1],axis[i])

class Binding(Predicate):
	"""docstring for Binding"""
	def __init__(self, obj,subj,direct=True,evidenceSentence=None):
		super(Binding, self).__init__("Binding",obj,subj,direct,evidenceSentence)

		
cmd=generate_variable_command("""ifng,ifngr,stat1,socs1,il-4,il-4r,stat6,gata3,t-bet,il-12,il-12r,stat4,il-18,il-18r,irak,ifnb,ifnbr,il21,stat3,tgfb,
il-1b, il-12, il-4, il-6, il-5, il-18, il-2, il-21, il-23,il-17,il-23r
il-2r
smad3, nfat,foxp3,stat5,jak3,foxp3prom
il17a,il17f,rorgt
""")
exec(compile(cmd,"","exec"))

Binding(stat3,il21,True,8)
Activation(il_6,stat3)
Activation([smad3,nfat],foxp3)
Activation(tgfb,smad3) #very old result
SignalingAxis([il_2,il_2r,jak3,stat5])
Binding(stat5,il17a)
Binding(stat3,il17a)
Binding(stat5,foxp3prom)
