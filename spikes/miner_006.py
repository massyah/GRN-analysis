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

from moduleTest import *

tc = ColorANSI.TermColors()



pubmed_files=[f for f in os.listdir(".") if f.startswith("entrez ")]
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
			t=Term(tt[0],alternatives)
		else:
			continue
	#parse pubmed results

	for pfile in pubmed_files:
		print "parsing",pfile
		txt=open(pfile,"r")
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


def omg_graph():
	existing_edges=[]
	for p in allPredicates:

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
	for t in [il_21,il_6,il_27,il_2,tcr,il_23,tgfb]:
		t.omg.rank_group.set(-2)
	for t in [il21,foxp3prom,il21,il17a,il17f,il21r,ifngprom]:
		t.omg.rank_group.set(-1)
	for t in [stat1,stat3,stat4,stat5]:
		t.omg.rank_group.set(0)
	for t in [il_2r,il_21r,il_23r,il_6r,il_1br]: #for receptors
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
	

global allPublications, allSentences, allTerms, allTermsRe, allPredicates, uid, sentencesUid, evidencesUid, pubmedIdTopub, nxG,cmd
print "default values"
allPublications={}
allSentences=[]
allTerms={}
allTermsRe=None
allPredicates=[]
uid=0
sentencesUid=0
evidencesUid=0
pubmedIdTopub={}
nxG=None	
parse_file()

cmd=generate_variable_command("""ifng,ifngr,stat1,socs1,il-4,il-4r,stat6,gata3,t-bet,il-12,il-12r,stat4,il-18,il-18r,irak,ifnb,ifnbr,il21,stat3,tgfb,ahr,
il-1b, il-12, il-4, il-6, il-5, il-18, il-2, il-21, il-23,il-17,il-23r,il-27
il-2r,il-6r,il-21r,il-27r
smad3, nfat,foxp3,stat5,jak3,foxp3prom,runx1,socs3,jak,stat1,jak1,jak2
il17a,il17f,rorgt,il21,il21r,il-1b,il-1br,il22
cd4p,tcr,nfkb,ap1,ets1,ifngprom,rora
myd88,
traf6,
jnk,
ikk,
jun,
fos,
jak2-jak1,
tyk2-jak2,
jak1-jak3,
ibp,irf4,smad4

""")
exec(compile(cmd,"","exec"))

mine=Publication("My findings","hayssam",2010)
ra=allTerms["_ra"]
il17_secretion=allTerms["il-17 secretion"]
il21_secretion=allTerms["il-21 secretion"]

oshea2009=allPublications["20090615 o'shea, john j"]

#signaling axis from the ILR
SignalingAxis([il_2,il_2r,jak1_jak3,stat5],oshea2009)
SignalingAxis([il_2,il_2r,jak1_jak3,stat3],oshea2009)
SignalingAxis([il_2,il_2r,jak1_jak3,stat1],oshea2009)
SignalingAxis([il_6,il_6r,jak2_jak1,stat3],7)
SignalingAxis([il_21,il_21r,jak1_jak3,stat3],7)
SignalingAxis([il_21,il_21r,jak1_jak3,stat1],7)
SignalingAxis([il_21,il_21r,jak1_jak3,stat5],7)
SignalingAxis([il_23,il_23r,tyk2_jak2,stat3],7)
SignalingAxis([il_23,il_23r,tyk2_jak2,stat4],14)
SignalingAxis([il_1b,il_1br,myd88,irak,traf6],27)



Activates([smad3,nfat],foxp3,[43,18157133])
Activates(tgfb,smad3) #very old result
Binds(stat3,il17a)
Binds(stat3,il17f,8)
Binds(stat5,foxp3prom)
Activates(tgfb,foxp3) #from PMID 14676299
Upregulates(t_bet,il_23r,1)
Activates(tcr,nfat,2)
Binds(nfat,il17a,3)
Activates(tcr,nfkb,4)
Activates(tcr,ap1,4)
Complexes(ap1,nfat,5)
Binds(stat4,ifng)
Binds(stat3,il21,9)
Upregulates(stat3,il_23r,10)

Downregulates(socs3,stat3,11) #also see <38>, socs3 negreg il-23-mediated stat3 signaling, does it holds for other il-n mediated stat3 activation?

Activates(ifng,stat1,12)
Upregulates([il_27,ifng],t_bet,13)
d=BindsRepress(stat5,il17a,15) #seem validated in pmid 17363300
d.hypothesis=True

Activates(allTerms["il-2 secretion"],il_2)
Inhibits(ets1,allTerms["il-2 secretion"],19)

Binds(stat5,foxp3prom)
Inhibits(foxp3,rorgt)
Inhibits(rorgt,foxp3)
Upregulates(il21r,il_21r)
Upregulates(foxp3prom,foxp3)
Activates(stat3,il21r,16)

Downregulates(stat5,rorgt,[17,20696842])

Binds([ets1,t_bet],ifngprom,18)
Upregulates(ifngprom,ifng,direct=True)

Activates(tgfb,smad3,20)
Activates([tgfb,il_6],rorgt,21)
Activates([tgfb,il_2],foxp3,21)

Upregulates([tgfb,il_6],il_23r,22)
Inhibits(stat3,foxp3,23)
Downregulates(tgfb,il_23r,24)
Downregulates(tgfb,t_bet,25)
Downregulates(tgfb,allTerms["il-2 secretion"],25)
Activates([il_1b,il_6],rorgt)
Activates(traf6,jnk,28)
Activates(traf6,ikk,28)
Phosphorylates(jnk,jun,29)
Activates([fos,jun],ap1,29)
Activates(ikk,nfkb,30)
Binds(rorgt,il17a,31)
Upregulates(stat3,rorgt,32)
Activates([stat3,rorgt],il17_secretion,33)
Activates(il17a,il17_secretion)
Activates(il17f,il17_secretion)
RequiredToActivate([rora,rorgt],il17_secretion,34)
Downregulates(ra,rorgt,35)
Upregulates(ra,foxp3,35)
BindsRepress(foxp3,il17a,36)
Downregulates(foxp3,il17_secretion,37)
Downregulates(foxp3,il21_secretion,37)
Upregulates(il21,il21_secretion)
Sequestrates(ibp,irf4,39)
Upregulates([irf4,stat3],rorgt,40)
Binds([runx1,rorgt],il17a,41)
Sequestrates(foxp3,runx1,17377532) #yields runx1 that can't dimerize with rorgt
Activates(tcr,ibp,15249594)
Upregulates([il_21,rorgt,stat3],il_23r,44)
Inhibits(smad3,il_2,1388)
build_nx_graph()
#according to thieffry, NFAT is required for almost everything. e.g. for foxp3

