#!/usr/bin/env python
# encoding: utf-8

#spike document how REgepx can be used to get exceprts

from Foundation import *
from AppKit import *
from Cocoa import *
from Quartz import *

import unicodedata
import os,math,re,time
from appscript import *

skim=app("Skim")

# documents={}
# page_seps={}
# from IPython.Debugger import Tracer; debug_here = Tracer()

from IPython import ColorANSI
from IPython.genutils import Term
tc = ColorANSI.TermColors()

#dehyphenation
hyphenation_re=re.compile("\w+-\s\w+")
termsWithHyphens=["IL-","TGF-β","STAT-","JAK-","IFN-","IRF\d*-","JAK\d","TRAF-"]
termsWithHyphens_re=re.compile("\w+|".join(termsWithHyphens)+"\w+")

def dehyphenated_version(match):
	w=match.group(0)
	test=w.replace("- ","-")
	if termsWithHyphens_re.match(test):
		res=test
	else:
		res=w.replace("- ","")
	return res
	
def dehyphenate(pdfText):
	return hyphenation_re.sub(dehyphenated_version,pdfText)

# allTerms={
# 
# }


def highlight(self,iPython=True):
	outputString=copy.copy(self.string)
	if iPython:
		return str(self.uid)+":\t"+allTermsRe.sub(pyHighlighter_call_back,outputString)
	else:
		return str(self.uid)+":\t"+allTermsRe.sub(highlighter_call_back,outputString)

def pyHighlighter_call_back(matchobj):
	return tc.Red+matchobj.group(0)+tc.Normal

def add_document(path,title,bodyText):
	documents[path]=(title,bodyText)
	print "added",len(bodyText),"chars doc"


search_results=[]
def ft_search(terms,expansion=False):
	global search_results
	print "Searching with",type(terms),terms
	if type(terms)!=list:
		terms=parse_query(terms,expansion)
	if len(terms)<1:
		print "query too short"
		return [],terms
	terms_re=build_re(terms)

	r=[]
	search_results=[]

	for path,doc in documents.items():
		bodyText=doc[1].lower()
		if terms_re.search(bodyText):
			r.append((path,doc))

	for path,doc in r:
		# The text argument to highlight is the stored text of the title
		text = doc[1]
		for ex in highlight_extracts(path,text,terms):
			search_results.append((path,ex[0],ex[1],ex[2],ex[3],ex[4],ex[5]))
	search_results.sort(key=lambda x:x[1],reverse=True)
	for i in range(min(40,len(search_results))):
		r=search_results[i]
		print i,":",r[0],r[2],"\n",r[4]
		print "------"*7
	# return search_results[:100],terms

hits=[]
hitScore=[]
def highlight_and_register_hit(m):
	global hits
	hits.append((m.span()[0],m.group(0)))
	return tc.Red+m.group(0)+tc.Normal

def highlight_extracts(path,docText,terms):
	global hits,hitScore
	hits=[]
	search_re=build_re(terms)
	highlited=search_re.sub(highlight_and_register_hit,docText)
	#we score the hits
	hitScore=[]
	for i in range(len(hits)):
		score=0
		#max span is 300
		differentKw=set()
		j=i
		while j<len(hits) and hits[j][0]-hits[i][0]<=300:
			score+=1
			differentKw.add(hits[j][1].lower())
			j+=1
		# print differentKw,score*math.exp(len(differentKw))
		hitScore.append((hits[i][0],score*math.exp(len(differentKw)**2)))
	#get the 30 best hits
	hitScore.sort(key=lambda x:x[1],reverse=True)
	best=hitScore[:30]
	#merge neighbors hits
	best.sort(key=lambda x:x[0])
	merged=[]
	i=0
	while i<len(best):
		j=i+1
		while j<len(best) and (best[j][0]-best[j-1][0]<=300):
			j+=1
		merged.append((best[i][1],best[i][0],best[j-1][0]+300))
		i=j
	merged.sort(key=lambda x:x[0],reverse=True)
	res=[]
	if path in page_seps:
		pageSeps=page_seps[path]
	else:
		pageSeps=None
	for h in merged:
		start=max(0,h[1]-60)
		excerpt=docText[start:h[2]+300]
		excerpt=search_re.sub(lambda x:tc.Red+x.group(0)+tc.Normal,excerpt)
		i=0
		if pageSeps:
			while i<len(pageSeps) and pageSeps[i]<h[1]:
				i+=1 
			if i!=len(pageSeps):
				i+=1
		res.append((h[0],i,h[1],excerpt,start,h[2]+300))
	return res


def get_page_seps(text):
	seps=[]
	idx=0
	while idx!=-1:
		idx=text.find("\n",idx+1)
		seps.append(idx)
	return seps

def add_pdf_doc(title,pdfPath):
	global page_seps
	pdfDoc=PDFDocument.new()
	pdfDoc.initWithURL_(NSURL.fileURLWithPath_(pdfPath))
	pdfText=dehyphenate(pdfDoc.string())
	#we update the page sep mapping based on linefeed chars
	page_seps[pdfPath]=get_page_seps(pdfText)
	add_document(pdfPath,title,pdfText)
	

def add_dt_doc(d):
	ext=os.path.splitext(d.path())[1]
	if ext in [".pdf"]:
		add_pdf_doc(d.name(),d.path())
	else:
		add_document(u"1",d.path(),d.name(),d.plain_text())
	print "Indexed ",d.name().encode("utf-8")

def parse_query(q,expansion=False):
	# print "Parsing query",q
	query_re=re.compile("(\".+?\")|([\w\-_]+)",re.UNICODE)
	terms=[]
	for m in query_re.finditer(q):
		for res in m.groups(1):
			if res != None and type(res)!=int:
				if res[0]=="\"":
					res=res[1:-1]
				if expansion and res in allTerms:
					res=allTerms[res]
					terms.append(res.name)
					terms.extend(res.alternativeNames)
				else:
					terms.append(res)
	return terms

def open_result(idx):
	res=search_results[idx]
	if res[0].isdigit():
		#a pubmed doc
		os.system("open http://www.ncbi.nlm.nih.gov/pubmed/%s"%(res[0]))
	else:
		os.system("open \"%s\""%res[0])
		time.sleep(0.5)
		skim.documents[1].go(to=skim.documents[1].pages[res[2]])
		#select
		skim.documents[1].selection.set([skim.documents[1].text.characters[res[5]:res[6]]])
	

def index_docs():
	for d in docs:
		add_pdf_doc(d,d)
		
def index_abstracts():
	for p in allPublications.values():
		add_document(p.pubmedID + " "+p.date + " "+p.title,p.title,p.title+", "+p.date+"\n"+p.firstAuthor+"\n"+p.abstract)

def search_open(q):
	ft_search(q,True)
	inp=""
	while True:
		inp=raw_input().strip()
		if len(inp)<1:
			return
		try:
			idx=int(inp)
		except:
			print "error"
			continue
		open_result(idx)
	
		
docs=[
"/Volumes/Gauss/Documents/Library/Human th17/32.pdf",
"/Volumes/Gauss/Documents/Library/Human th17/207_ftp.pdf",
"/Volumes/Gauss/Documents/Library/Human th17/299.full.pdf",
"/Volumes/Gauss/Documents/Library/Human th17/469_ftp.pdf",
"/Volumes/Gauss/Documents/Library/Human th17/537.full.pdf",
"/Volumes/Gauss/Documents/Library/Human th17/637_ftp.pdf",
"/Volumes/Gauss/Documents/Library/Human th17/661_ftp.pdf",
"/Volumes/Gauss/Documents/Library/Human th17/959.full.pdf",
"/Volumes/Gauss/Documents/Library/Human th17/1291.full.pdf",
"/Volumes/Gauss/Documents/Library/Human th17/1742.pdf",
"/Volumes/Gauss/Documents/Library/Human th17/2213.pdf",
"/Volumes/Gauss/Documents/Library/Human th17/2983.pdf",
"/Volumes/Gauss/Documents/Library/Human th17/3610.pdf",
"/Volumes/Gauss/Documents/Library/Human th17/9741.pdf",
"/Volumes/Gauss/Documents/Library/Human th17/17034.full.pdf",
"/Volumes/Gauss/Documents/Library/Human th17/22866_ftp.pdf",
"/Volumes/Gauss/Documents/Library/Human th17/7100114a.pdf",
"/Volumes/Gauss/Documents/Library/Human th17/annurev.immunol.25.022106.141711.pdf",
"/Volumes/Gauss/Documents/Library/Human th17/annurev.immunol.25.022106.141557.pdf",
"/Volumes/Gauss/Documents/Library/Human th17/ar2392.pdf",
"/Volumes/Gauss/Documents/Library/Human th17/fulltext.pdf",
"/Volumes/Gauss/Documents/Library/Human th17/Gut-2009-Rovedatti-1629-36.pdf",
"/Volumes/Gauss/Documents/Library/Human th17/Int. Immunol.-2008-Annunziato-1361-8.pdf",
"/Volumes/Gauss/Documents/Library/Human th17/Int. Immunol.-2009-Minegishi-105-12.pdf",
"/Volumes/Gauss/Documents/Library/Human th17/Int. Immunol.-2010-Ashino-503-13.pdf",
"/Volumes/Gauss/Documents/Library/Human th17/J Biochem-2010-Yoshimura-781-92.pdf",
"/Volumes/Gauss/Documents/Library/Human th17/J Exp Med-2007-Annunziato-1849-61.pdf",
"/Volumes/Gauss/Documents/Library/Human th17/J Exp Med-2008-Cosmi-1903-16.pdf",
"/Volumes/Gauss/Documents/Library/Human th17/J Exp Med-2008-de Beaucoudrey-1543-50.pdf",
"/Volumes/Gauss/Documents/Library/Human th17/J Exp Med-2008-Rao-3145-58.pdf",
"/Volumes/Gauss/Documents/Library/Human th17/J Immunol-2003-Sundrud-3542-9.pdf",
"/Volumes/Gauss/Documents/Library/Human th17/J Immunol-2007-Peluso-732-9.pdf",
"/Volumes/Gauss/Documents/Library/Human th17/J Immunol-2007-Sato-7525-9.pdf",
"/Volumes/Gauss/Documents/Library/Human th17/J Immunol-2008-Lim-122-9.pdf",
"/Volumes/Gauss/Documents/Library/Human th17/J Immunol-2008-Pène-7423-30.pdf",
"/Volumes/Gauss/Documents/Library/Human th17/J Immunol-2009-Loures-1279-90.pdf",
"/Volumes/Gauss/Documents/Library/Human th17/J Immunol-2009-Nichols-7119-30.pdf",
"/Volumes/Gauss/Documents/Library/Human th17/J Immunol-2009-Tangye-21-8.pdf",
"/Volumes/Gauss/Documents/Library/Human th17/J Immunol-2009-Wang-4119-26.pdf",
"/Volumes/Gauss/Documents/Library/Human th17/J Immunol-2010-Crome-3199-208.pdf",
"/Volumes/Gauss/Documents/Library/Human th17/J Immunol-2010-Sellge-2076-85.pdf",
"/Volumes/Gauss/Documents/Library/Human th17/J. Biol. Chem.-2004-Liu-52762-71.pdf",
"/Volumes/Gauss/Documents/Library/Human th17/j.1600-065X.2008.00628.x.pdf",
"/Volumes/Gauss/Documents/Library/Human th17/j.1600-065X.2008.00703.x.pdf",
"/Volumes/Gauss/Documents/Library/Human th17/j.1600-065X.2008.00705.x.pdf",
"/Volumes/Gauss/Documents/Library/Human th17/j.1600-065X.2008.00712.x.pdf",
"/Volumes/Gauss/Documents/Library/Human th17/nature06878.pdf",
"/Volumes/Gauss/Documents/Library/Human th17/nature07021.pdf",
"/Volumes/Gauss/Documents/Library/Human th17/ni.1610.pdf",
"/Volumes/Gauss/Documents/Library/Human th17/ni.1613.pdf",
"/Volumes/Gauss/Documents/Library/Human th17/ni.1767.pdf",
"/Volumes/Gauss/Documents/Library/Human th17/ni0608-588.pdf",
"/Volumes/Gauss/Documents/Library/Human th17/ni0907-903.pdf",
"/Volumes/Gauss/Documents/Library/Human th17/ni1449.pdf",
"/Volumes/Gauss/Documents/Library/Human th17/ni1460.pdf",
"/Volumes/Gauss/Documents/Library/Human th17/ni1467.pdf",
"/Volumes/Gauss/Documents/Library/Human th17/ni1488.pdf",
"/Volumes/Gauss/Documents/Library/Human th17/ni1496.pdf",
"/Volumes/Gauss/Documents/Library/Human th17/ni1497.pdf",
"/Volumes/Gauss/Documents/Library/Human th17/ni1500.pdf",
"/Volumes/Gauss/Documents/Library/Human th17/nm0201_245.pdf",
"/Volumes/Gauss/Documents/Library/Human th17/nm1585.pdf",
"/Volumes/Gauss/Documents/Library/Human th17/nm1651.pdf",
"/Volumes/Gauss/Documents/Library/Human th17/nri1688.pdf",
"/Volumes/Gauss/Documents/Library/Human th17/nri2295.pdf",
"/Volumes/Gauss/Documents/Library/Human th17/PNAS-2008-Miyahara-15505-10.pdf",
"/Volumes/Gauss/Documents/Library/Human th17/PNAS-2009-Mars-6238-43.pdf",
"/Volumes/Gauss/Documents/Library/Human th17/PNAS-2009-Voo-4793-8.pdf",
"/Volumes/Gauss/Documents/Library/Human th17/Rheumatology-2009-Nistala-602-6.pdf",
"/Volumes/Gauss/Documents/Library/Human th17/science-1.pdf",
"/Volumes/Gauss/Documents/Library/Human th17/science-2.pdf",
"/Volumes/Gauss/Documents/Library/Human th17/science-3.pdf",
"/Volumes/Gauss/Documents/Library/Human th17/science-4.pdf",
"/Volumes/Gauss/Documents/Library/Human th17/science-5.pdf",
"/Volumes/Gauss/Documents/Library/Human th17/Science-2007-Mucida-256-60.pdf",
"/Volumes/Gauss/Documents/Library/Human th17/science.pdf",
"/Volumes/Gauss/Documents/Library/Human th17/sdarticle-1.pdf",
"/Volumes/Gauss/Documents/Library/Human th17/sdarticle.pdf",
"/Volumes/Gauss/Documents/Library/Human th17/Unknown-1.pdf",
"/Volumes/Gauss/Documents/Library/Human th17/Unknown.pdf"]


textTest=[
"""Enhanced glycolysis is important for oncogenesis and for the survival and proliferation of cancer cells in the tumor microenvironment. Recent studies have also shown that proinflammatory cytokine signaling, such as that mediated by nuclear factor kappaB and signal transducer and activator of transcription 3 (STAT3), is important for the generation of inflammation-associated tumors. However, the link between inflammation and enhanced glycolysis has not been identified. In the present study, we found that the proinflammatory cytokine interleukin (IL)-6 enhanced glycolysis in mouse embryonic fibroblasts and human cell lines. Moreover, STAT3 activated by IL-6 enhanced the expression of the glycolytic enzymes hexokinase 2 and 6-phosphofructo-2-kinase/fructose-2,6-bisphosphatase-3 (PFKFB3). Ectopic expression of PFKFB3 enhanced glycolysis, suggesting that the IL-6-STAT3 pathway enhances glycolysis through the induction of these enzymes. Our findings may provide a novel mechanism for inflammation-associated oncogenesis.""",
"""d to be mutually exclusive. Importantly, in the presence of IL-6, TGF-β induced development of Treg cells was blocked, whereas blockade of IL-6 permitted development of Foxp3+ Treg cells, suggesting that IL-6 inhibited Treg development while enhancing Th17 development induced by TGF-β [17,24,25]. In a recent study, however, an IL-6 independentProperties of murine Th17 cells and their role in protection and immunopathology The master regulator that directs the differentiation program of Th17 cells is the orphan retinoid nuclear receptor (ROR)γt, whereas neither GATA-3 nor T-bet are required for this function [29]. More recently, it was found that Th17 cells express high levels of anThe major functions of cytokines produced by Th17 cells is to chemIL-6, which is dependent on STAT-3 [30].""",
"""n of FoxP3+ Treg induced by TGF-b. Furthermore, analysis of IL-6-deficient mice revealed that the ability of IL-21 to generate Th17 cells from CD4+ T cells was independent of IL-6. Thus, IL-21, like IL-6, can dictate the generation of Th17 versus Treg. For both IL-21 and IL-6, this switch seems to be mediated by STAT3 and RORgt.6–8 These papers also demonstrated that IL-6 or IL-21 could induce Th17 cells themselves to produce IL-21. Such endogenous production of IL-21 by Th17 cells appeared to be biologically significant, because the number of IL-17-producing cells generated by TGF-b and IL-6 was reduced in the absence of IL-21/IL-21R signalling.6,7 Thus, IL-6 can elicit IL-21 production by CD4+ T cells, which then functions in an autocrine loop to amplify the Th17 response in a similar way to IL-4 for Th2 and IFN-g for Th1 cells (Figure 1). Furthermore, both IL-21 and IL-6 upregulated IL-23R, thereby primiTargeting molecules involved in the generation, maintenance and effector function of Th17 cells, such as TGF-b, IL-6, IL-23 and IL-17, have been proposed as novel therapeutics for human inflammatory disorders. These findings on the role of IL-21 in the biology of Th17 cells6–8 add IL-21 to this list of potential targets. In particular, since IL-21 functions in a self-amplifying loop (Figure 1) on established Th17 cells, IL-21 may be a better target for inhibiting the pathogenic inflammatory response than IL-6. Recent reports showing that neutralizing IL-21 in murine models of lupus or rheumatoid arthritis ameliorated disease severity17,18 demonstrate the potential utility of IL-21targeted therapies for """,
"""IL-22-producing T cells (%) IL-17-producing T cells (%) IL-22 IL-22 IFN-γ-producing T cells (%) 100101 102103 104 100101 102103 104 100101 102103 104 103Medium 4 IL-6 10IL-6 + TNF 100 10010 1 0 1 1010 2 2 1010 3 10 104 IL-6 + TNF + IL-1βIL-6 + TNF + TGF-β (5) IL-6 + TNF + TGF-β (0.5)IL-6 + TNF + TGF-β (0.05) 10110 100101 10210 2 10310 4 3104 MediumIL-6 IL-6 + TNFIL-6 + TNF + IL-1β IL-6 + TNF + TGF-β (5)IL-6 + TNF + TGF-β (0.5) IL-6 + TNF + TGF-β (0.05) 100101 102103 104 MediumIL-6 IL-6 + TNFIL-6 + TNF + IL-1β IL-6 + TNF + TGF-β (5)IL-6 + TNF + TGF-β (0.5) IL-6 + TNF + TGF-β (0.05) a b 40302010 Medium IL-68.1 0.2 0.8 5.9 2.411.0 IL-6 + TNF10.0 3.0 0.4 IL-6 + TNF + IL-1β5.8 9.3103 0.9 102101100 104 102101100 000 104103102101100 IL-17104103102101100 IFN-γ P = NSP = 0.004 6.1 0.40.2 2.9 3.713.2 104103102101100 104103102101100 104103102101100 104103102101100 104 10.0 3.04.1 103 8.5 6.65.5 40 30 20 10 P = 0.003 P = NS 60 40 20 P = NSP = NS Figure 4 IL-6 and TNF promote TH-22 differentiation. (a) Cytokine production by T naive CD4+CD45RA+CCR7+CD45RO–CD25– T cells stimulated in vitro for 5 d with beads coated with anti-CD3and anti-CD28 in the absence (Medium) or presence of various combinations of IL-6, TNF and ILb(above plots) and expanded with IL-2 for an additional 7 d, followed by intracellular cytokine stainingafter restimulation for 5 h with PMA and ionomycin. Data are representative of eight independentexperiments. (b) Frequency of T cells producing IL-22, IL-17 or IFN-g among naive T cells primed withvarious combinations of IL-6, TNF and IL-1b, plus TGF-b at various concentrations (horizontal axis, inng/ml). Each symbol represents one donor; small horizontal bars indicate the mean. P values, Student’st-test. Data are from five to eight independent experiments."""
]

