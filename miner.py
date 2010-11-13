#!/usr/bin/env python
# encoding: utf-8
import re,copy,os
import omg_interface as omg
import pubmed_get as pm
import AppKit
import networkx as nx

from IPython import ColorANSI
from IPython.genutils import Term
from Bio import Entrez,Medline

tc = ColorANSI.TermColors()



# pubmed_files=["entrez 19379825.txt","entrez 18368049.txt"]

evid_re=re.compile("<([0-9]+)>")

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
		self.pubmedFileCollection=None
		self.fulltextFile=None
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
		self.evidenceUid=None
		self.toUnderstand=False
		sentencesUid+=1
		allSentences.append(self)
	def term_recognition(self):
		global evidencesUid
		if self.string.startswith("<"):
			if evid_re.search(self.string):
				self.evidenceUid=int(evid_re.findall(self.string)[0])
				evidencesUid=max(self.evidenceUid+1,evidencesUid)
			else:
				self.toUnderstand=True
		for m in allTermsRe.findall(self.string):
			if m.lower() in allTerms:
				self.terms.append(allTerms[m.lower()])
			elif "_"+m.lower() in allTerms:
				self.terms.append(allTerms["_"+m.lower()])
			else:
				self.terms.append(allTerms["_"+m.lower()])
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


def parse_file(aTermFile=None):
	if aTermFile!=None:
		termFile=aTermFile
	lines=open(termFile).readlines()
	print "#",len(lines),"lines in term file"
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
			t=Term(tt[0],alternatives)
		else:
			continue
	#parse pubmed results

	for pfile in pubmed_files:
		print "parsing",pfile
		txt=open(project+"/"+pfile,"r")
		records=Medline.parse(txt)
		for r in records:
			if "AB" not in r:
				continue
			if "FAU" in r:
				firstAuthor=r["FAU"][0]
			elif "AU" in r:
				firstAuthor=r["AU"][0]
			else:
				firstAuthor="Unknown"
			date=r["DA"]
			title=r["TI"]
			pmid=r["PMID"]
			if pmid not in pubmedIdTopub:
				pub=Publication(title,firstAuthor,date,pmid)
				pub.pubmedFileCollection=pfile
				pubmedIdTopub[pmid]=pub
				for li in r['AB'].split(". "):
					Sentence(pub,li)
				fullTxtFile="pmid %s.txt"%(pmid)
				if os.path.isfile(fullTxtFile):
					print "should parse full text ",fullTxtFile
					pub.fulltextFile=fullTxtFile
					fullTxt=open(fullTxtFile).read()
					#clean 
					fullTxt=fullTxt.replace("Fig .","Fig ")
					fullTxt=fullTxt.replace("et al.","et al ")
					fullTxt=fullTxt.replace("e.g.","eg ")
					fullTxt=fullTxt.replace("vs.","vs")
					fullTxt=fullTxt.replace("ref.","ref ")
					fullTxt=unicode(fullTxt,"utf-8")
					for li in fullTxt.split("."):
						Sentence(pub,li.strip())

	for s in allSentences:
		s.term_recognition()
	
def generate_variable_command(str,sep=","):
	cmd=""
	for l in str.split("\n"):
		if len(l)<2:
			continue
		for v in l.split(sep):
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
	
def print_graph():
	for p in allPredicates:
		print p.__str__()

def annotate_sentence_with_evidence_uid(tgtUid):
	global evidencesUid
	#search for sentence in the file according ot its uid
	for s in allSentences:
		if s.uid==tgtUid:
			tgtSent=s
			break

	os.system("open -a TextMate \"%s\""%(tgtSent.publication.fulltextFile))

	#put the sentence on the clip
	p=AppKit.NSPasteboard.pasteboardWithName_(AppKit.NSFindPboard)
	p.declareTypes_owner_([AppKit.NSStringPboardType],None)
	p.setString_forType_(tgtSent.string,AppKit.NSStringPboardType)
	#generate the evidence uid and put it on the pboard
	uidStr="<%d>"%(evidencesUid)
	p=AppKit.NSPasteboard.generalPasteboard()
	p.declareTypes_owner_([AppKit.NSStringPboardType],None)
	p.setString_forType_(uidStr,AppKit.NSStringPboardType)
	evidencesUid+=1

class Predicate(object):
	"""docstring for Predicate"""
	def __init__(self, name,obj,subj,evidenceSentence=None,direct=False):
		global allTerms, uid
		super(Predicate, self).__init__()
		self.name = name
		self.obj = obj
		self.subj = subj
		self.uid=uid
		self.direct=direct
		self.hypothesis=False
		self.evidenceSentence=evidenceSentence
		uid+=1
		allPredicates.append(self)
	def __str__(self):
		return "%s -> %s"%(self.obj.name,self.subj.name)





#first mappings
class Inactivates(Predicate):
	def __init__(self, obj,subj,evidenceSentence=None,direct=False):
		super(Inactivates, self).__init__("Inactivates",obj,subj,evidenceSentence,direct)

class Activates(Predicate):
	def __init__(self, obj,subj,evidenceSentence=None,direct=False):
		super(Activates, self).__init__("Activates",obj,subj,evidenceSentence,direct)
	def __str__(self):
		return "%s -+> %s"%(self.obj.name,self.subj.name)

class Inhibits(Predicate):
	def __init__(self, obj,subj,evidenceSentence=None,direct=False):
		super(Inhibits, self).__init__("Inhibits",obj,subj,evidenceSentence,direct)
	def __str__(self):
		return "%s -| %s"%(self.obj.name,self.subj.name)
	
def SignalingAxis(axis,evidenceSentence=None):
	for i in range(1,len(axis)):
		Activates(axis[i-1],axis[i],evidenceSentence)

class Binds(Predicate):
	def __init__(self, obj,subj,evidenceSentence=None,direct=True):
		super(Binds, self).__init__("Binds",obj,subj,evidenceSentence,direct)

class BindsRepress(Binds):
	def __init__(self, obj,subj,evidenceSentence=None,direct=True):
		super(BindsRepress, self).__init__(obj,subj,evidenceSentence,direct)

class Downregulates(Predicate):
	def __init__(self, obj,subj,evidenceSentence=None,direct=False):
		super(Downregulates, self).__init__("Downregulates",obj,subj,evidenceSentence,direct)

class Upregulates(Predicate):
	def __init__(self, obj,subj,evidenceSentence=None,direct=False):
		super(Upregulates, self).__init__("Upregulates",obj,subj,evidenceSentence,direct)
	
class Complexes(Predicate):
	def __init__(self, obj,subj,evidenceSentence=None,direct=True):
		super(Complexes, self).__init__("Complexes",obj,subj,evidenceSentence,direct)
		
class Phosphorylates(Predicate):
	def __init__(self, obj,subj,evidenceSentence=None,direct=True):
		super(Phosphorylates, self).__init__("Phosphorylates", obj,subj,evidenceSentence,direct)

class RequiredToActivate(Activates):
	def __init__(self, obj,subj,evidenceSentence=None,direct=False):
		super(RequiredToActivate, self).__init__(obj,subj,evidenceSentence,direct)

class Sequestrates(Inactivates):
	def __init__(self, obj,subj,evidenceSentence=None,direct=False):
		super(Sequestrates, self).__init__(obj,subj,evidenceSentence,direct)


def omg_graph(predicates=None):
	global ranks
	if predicates==None:
		predicates=allPredicates
	existing_edges=[]
	for p in predicates:

		if type(p) in [Activates]:
			f=omg.add_activation
		elif type(p) in [Upregulates]:
			f=omg.add_upregulation
		elif type(p) in [Downregulates,Inhibits,Inactivates,Sequestrates]:
			f=omg.add_inhibition
		elif type(p) in [BindsRepress]:
			f=omg.add_inhibition
		elif type(p) in [Binds]:
			f=omg.add_binding
		else:
			f=omg.add_activation
			

		if type(p.obj)==list:
			a=omg.add_and_node(p.obj)
			f(a,p.subj)
		else:
			if (p.obj,p.subj) in existing_edges:
				continue
			existing_edges.append((p.obj,p.subj))
			f(p.obj,p.subj)
	#assign rank groups
	for t in ranks["sign"]:
		if t.omg != None:
			t.omg.rank_group.set(-2)
	for t in ranks["proms"]:
		if t.omg != None:
			t.omg.rank_group.set(-1)
	for t in ranks["stats"]:
		if t.omg != None:	
			t.omg.rank_group.set(0)
	for t in ranks["receptors"]: #for receptors
		if t.omg != None:
			t.omg.rank_group.set(-50)
	omg.layout()

def build_nx_graph(onlyPositive=True):
	global nxG
	uniqueAndId=0
	nxG=nx.DiGraph()
	for p in allPredicates:
		if (onlyPositive and type(p) in [Activates, Binds, Complexes, Phosphorylates, RequiredToActivate, Upregulates])\
		or not onlyPositive:
				if type(p.obj)==list:
					uniqueAndNodeStr="AND %d"%(uniqueAndId)
					uniqueAndId+=1
					for obj in p.obj:
						nxG.add_edge(obj.name,uniqueAndNodeStr,{"uid":p.uid})
					nxG.add_edge(uniqueAndNodeStr,p.subj.name,{"uid":p.uid})
				else: #not list
					nxG.add_edge(p.obj.name,p.subj.name,{"uid":p.uid})

def predicates_with_uid(uid):
	return [p for p in allPredicates if p.uid==uid]

def sentence_with_evidence_uid(uid):
	return [s for s in allSentences if s.evidenceUid==uid]

def shortest_path(src,tgt):
	p=nx.shortest_path(nxG,src.name,tgt.name)
	justif=[]
	print "+->+".join(p)
	#assoc with predicates
	for i in range(len(p)-1):
		print nxG[p[i]][p[i+1]]
		pred=predicates_with_uid(nxG[p[i]][p[i+1]]["uid"])[0]
		justif.append(sentence_with_evidence_uid(pred.evidenceSentence)[0].string)
	print "\n".join(justif)

def term_statistics_in_sentences(sents):
	termN={}
	for s in sents:
		for t in s.terms:
			if t.name not in termN:
				termN[t.name]=0
			termN[t.name]+=1
	termN=sorted(termN.iteritems(), key=lambda x:x[1], reverse=True)
	for k,v in termN:
		print k,v