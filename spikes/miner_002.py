#!/usr/bin/env python
# encoding: utf-8
import re,copy

from IPython import ColorANSI
from IPython.genutils import Term
tc = ColorANSI.TermColors()


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
def pyHighlighter_call_back(matchobj):
	return tc.Red+matchobj.group(0)+tc.Normal
	
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
	lines=open("articles.txt").readlines()
	print len(lines)
	lastPub=None
	global allPublications,allTerms,allSentences,allPredicates,uid
	uid=0
	allPublications={}
	allTerms={}
	allSentences=[]
	allPredicates=[]
	for l in lines:
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
			t=Term(tt[0],alternatives)
		elif l.startswith("#p"): #create a new pub
			l=l[2:].strip()
			lastPub=Publication(*l.split(","))
		elif l.startswith("/"):
			continue
		else:
			assert(lastPub)
			Sentence(lastPub,l)

			
	for s in allSentences:
		s.term_recognition()
	
def generate_variable_command(str):
	cmd=""
	vars=str.split(",")
	for v in vars:
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
print len(allSentences),len(allTerms),len(allPublications)

#some terms from mendoza

def print_graph():
	for p in allPredicates:
		print p.__str__()

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
		allPredicates.append(self)
	def __str__(self):
		return "%s -> %s"%(self.obj.name,self.subj.name)





#first mappings
class Activation(Predicate):
	def __init__(self, obj,subj):
		super(Activation, self).__init__("Activation",obj,subj)
	def __str__(self):
		return "%s -+> %s"%(self.obj.name,self.subj.name)

class Inhibition(Predicate):
	def __init__(self, obj,subj):
		super(Inhibition, self).__init__("Inhibition",obj,subj)
	def __str__(self):
		return "%s -| %s"%(self.obj.name,self.subj.name)
	
def SignalingAxis(axis):
	for i in range(1,len(axis)):
		Activation(axis[i-1],axis[i])

cmd=generate_variable_command("ifng,ifngr,stat1,socs1,il-4,il-4r,stat6,gata3,t-bet,il-12,il-12r,stat4,il-18,il-18r,irak,ifnb,ifnbr")
exec(compile(cmd,"","exec"))

SignalingAxis([ifng,ifngr,stat1,socs1]) #39,46,55,60
SignalingAxis([ifng,ifngr,stat1,t_bet]) #104
SignalingAxis([il_4,il_4r,stat6,gata3]) #114, ~123
SignalingAxis([il_12,il_12r,stat4]) #?
SignalingAxis([il_18,il_18r,irak])#?
SignalingAxis([ifnb,ifnbr,stat1]) #93
Activation(gata3,gata3) 
Inhibition(gata3,stat4) #126,127
Activation(t_bet,socs1) #?
Activation(t_bet,t_bet) #109
Inhibition(stat6,il_18r)#?
Inhibition(stat6,il_12r)#?
Inhibition(t_bet,gata3) #111
Activation(stat4,ifng) #?
Activation(irak,ifng) #?
Inhibition(socs1,ifngr) #~72, weak
