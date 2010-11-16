#!/usr/bin/env python
# encoding: utf-8

import os,itertools,unicodedata,sys

project="Th17 differentiation"
termFile=project+"/articles_th17.txt"
pubmed_files=[f for f in os.listdir("./"+project) if f.startswith("entrez ")]
pubmed_files=[]

model_name="th17.ginml"

from IPython.Debugger import Tracer; debug_here = Tracer()


import omg_interface as omg

global allPublications, allSentences, allTerms, allTermsRe, allPredicates, uid, sentencesUid, evidencesUid, pubmedIdTopub, nxG,cmd

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
all_edges_id=[]


termDefs=u"""
//Terms
//Cytokines
#t IL-17
#t IL-1β, IL-1B
#t IL-21
#t IL-23
#t IL-6
#t TGF-β, tgf beta, tgfb, TGF-beta,TGFbeta,TGFβ


//Receptors
#t IL-12 receptor β2, IL-12rβ2,IL-12R2,IL-12 receptor 2
#t IL-1βR, IL-1BR
#t IL-21R
#t IL-23R,IL23-r
#t IL-6R
#t TGF-βR, TGFBR

//SMADS
#t SMAD3
#t SMAD4, Co-SMAD

//Jaks & STATs
#t STAT-3,STAT3
#t JAK
#t JAK-3, JAK3
#t JAK-2, JAK2
#t JAK-1, JAK1

//Jak & Tyk dimers
#t Jak2-Jak1
#t Tyk2-Jak2
#t Jak1-Jak3

//Transcription factors
#t FOXP3
#t RUNX1, AML1, CBFA2, PEBP2A2, AMLCR1
#t RUNX2, AML2


#t T-bet, TBX21, TBET
//From http://www.genenames.org/data/hgnc_data.php?hgnc_id=11599

//Nuclear receptor transcription factors
#t ROR-γt,RORγt,rorgt, RORgammat, retinoic acid-related orphan receptor gammat
#t RORα,RORa

//Specific genes locii related to secretion of th17 cytokines
#t il17a,IL-17A, il-17a promoter
#t il17f, il-17f, IL-17F promoter
#t il6 
#t il21
#t il21r 
#t il23r
#t foxp3prom, FOXP3 promoter 
#t rorgtprom, RORγt promoter
"""

f=open("th17_model.terms","w")
f.write(termDefs.encode("utf-8"))
f.close()
parse_file("th17_model.terms")

cmd=generate_variable_command("""tgfb,foxp3,rorgt, jak2-jak1, tyk2-jak2, jak1-jak3,stat3,

il-1b,il-21, il-6,il-23,il-17 ,il-6
il-6r,il-21r,il-23r,il-1br,

il23r,il21,il21r,il17a,il17f
""")
exec(compile(cmd,"","exec"))

ranks={
"sign":[il_21,il_6,il_23,il_17,tgfb],
"stats":[stat3],
"proms":[il21,il21r,il17a,il17f],
"receptors":[il_21r,il_23r,il_6r,il_1br]
}




experimentalResults=[]
def domain(var):
	nonOne={tgfb:3,il23r:2,rorgt:2,stat3:2}
	if type(var) == unicode:
		var=allTerms[var.lower()]
	if var in nonOne:
		return nonOne[var]
	else:
		return 1

def domain_for_vars(*args):
	ranges=[range(domain(x)+1) for x in args]
	vals=[x for x in itertools.product(*ranges)]
	return vals


def input_xp_results(inpVar,outVar):
	#generate the input tuple values, 
	inpVar=map(lambda x:x.name,inpVar)
	outVar=map(lambda x:x.name,outVar)
	blankTable=""
	possible_input=domain_for_vars(*inpVar)
	print "\t".join(inpVar+outVar)
	for v in possible_input:
		blankTable+="\t".join(map(str,v))+"\t?"+"\n"
	return blankTable
	
all_results=[]
all_experiments=[]
class Experiment(object):
	"""docstring for Experiment"""
	def __init__(self, skimReference,inpVars,outVars,results):
		global all_experiments
		super(Experiment, self).__init__()
		self.skimReference = skimReference
		self.results = results
		self.relation=[]
		self.inpVars=inpVars
		self.outVars=outVars
		all_experiments.append(self)
		self.parse_results()
	def parse_results(self):
		global all_results
		if self.results==None:
			return
		for l in self.results.split("\n"):
			if len(l.strip())<1:
				continue
			if l[0].isalpha():
				continue
			res=tuple(zip(self.inpVars+self.outVars,l.split("\t")))
			inpDict,outDict=dict(res[:len(self.inpVars)]),dict(res[len(self.inpVars):])
			self.relation.append((inpDict,outDict))
			all_results.append((inpDict,outDict,self))
			
def known_results(inpVars,outVars):
	res=[]
	inpVars,outVars=set(inpVars),set(outVars)	
	for inp,out,exp in all_results:
		inpS,outS=set(inp.keys()),set(out.keys())
		if inpVars<=inpS and outVars<=outS:
			res.append((inp,out,exp))
	return res

def print_results(res):
	res.sort(key=lambda x:x[2])
	groups=itertools.groupby(res,lambda x:x[2])
	for k,g in groups:
		header=True
		for inp,out,exp in g:
			if header:
				print exp.skimReference
				print "\t".join(map(lambda x:x.name,inp.keys()+out.keys()))
				header=False

			print "\t".join(inp.values()+out.values())
		print "------"
		
		
#assumed interactions
SignalingAxis([il_6,il_6r,jak2_jak1,stat3])
SignalingAxis([il_21,il_21r,jak1_jak3,stat3])

# Activates(stat3, rorgt) #Should it be cooperative binding instead?


Inhibits(foxp3, rorgt,"12.5.1	Only the full-length isoform [of foxp3]... was found to co-precipitate with RORγt and to inhibit its function")


Activates(tgfb, foxp3)

Binds(rorgt, il17a,"4.6.2	RORγt has been found to bind the Il17A gene [81].")
Activates(tgfb,rorgt)

Binds(rorgt, il23r)

Binds(stat3,il23r,"4.3.3 Furthermore, IL-6, IL-21 and IL-23 can upregulate expression of the IL-23R, and this too appears to be Stat3-dependent [9,10,16].")
Binds(stat3, il21,"4.3.1 Stat3 appears to be a direct regulator of IL-21 production as it binds to the Il21 promoter [32].")


Binds(stat3,il17a, "4.3.2")
Binds(stat3,il17f,"4.3.2 The first evidence of the importance of Stat3 in regulating IL-17 production came from the demonstration that the Il17A-Il17F locus has multiple putative Stat binding sites. Using chromatin immunoprecipitation (ChIP) assays, it was determined that Stat3 directly binds to this promoter [48]")

Upregulates(il23r,il_23r)
Upregulates(il21,il_21)
Upregulates(il17a,il_17)
Upregulates(il17f,il_17)

Inhibits(stat3,foxp3,"4.5.1	Additionally, IL-6 represses Foxp3 in a Stat3-dependent manner [71,76].")

BindsRepress(foxp3,il17a,"""
4.6.6	Supporting this is evidence that Foxp3 binds to the Il17A promoter [80]. #No evidence in 80!!!
12.5.2	Furthermore, mutations that abolished Foxp3 DNA binding also blunted its ability to inhibit IL-17 expression, suggesting that Foxp3 inhibits expression of RORγt targets by directly binding to the transcription factor and by also acting upon the target gene (Liang Zhou et al., submitted).
""")

todos=[
"4.5.2	Also of note is that IL-17 and IL-22 appear to differ in their requirement for TGF-β [83].",
"4.5.3	[IL-1 posreg th17] This has been documented in both mouse and human cells where it synergizes with IL-6 and IL-23 [67,90].",
"4.5.4	Whereas IL-6 alone causes transient induction of Rorgt, the combination of IL-1 and IL-6 results in sustained induction.",
"4.6.1	Conversely, overexpression of RORγt promotes IL-17 expression; in this regard RORγt acts similarly to transcription factors such as T-bet and GATA3 and therefore has been proposed to be a ‘‘master regulator’’ for Th17 differentiation.",
"4.6.3	Overexpression of active Stat3 in Rorc-deficient cells resulted in poor IL-17 induction, indicating that Stat3 is necessary but not sufficient for IL-17 expression."
"4.6.4	However, overexpression of RORγt failed to induce IL-17 secretion, suggesting that these two transcription factors act in parallel to some extent [80].",
"4.6.5	Stat3 mediates IL-6 regulation of RORγt, but it has not been established that Stat3 binds the Rorc gene. #no ref?",
"4.6.7	IL-17 and IL-21 are down-regulated in Foxp3-expressing cells [104], but it is not known if Foxp3 also directly binds the Rorc or Il21 promoters.",
"diff between il17a and il17f?",
"learn wikipedia page of all the specie in the model, do not make a single mistake!"
"12.3.3	Like IL-21, IL-23R mRNA is induced after treatment of TCR-activated T cells with IL-6 or IL-21 and is also upregulated by forced expression of RORγt [26,29]. IL-23R expression was abrogated in IL-21R-deficient cells, suggesting that IL-21 is required as an intermediate in IL-6-mediated induction of IL-23R."
"12.5.3	In the presence of proinflammatory cytokines, TGF-- induced Foxp3 expression is greatly reduced and the RORγt level is further up-regulated, thereby favoring Th17 cell differentiation [10,26] #is this incompatible with an inactivation by STAT3? Is it more complex (on the expression level?)"
"1.2.3	In addition to the aforementioned transcription factors, interferon-regulatory factor 4 and T-bet are two relative newcomers to the scene, and these act as positive and negative regulators of Th17 commitment, respectively [51, 52]."
"5.7.1	TCR stimulation as well as exposure to IL-6 lead to downregulation and shedding of the IL-6Rα and thus reduce responsiveness to IL-6. TGF-β is necessary to maintain the responsiveness of T cells to IL-6."
"5.7.2	Deficiency of STAT3 in T cells abrogates the induction of Th17 cells and decreases the induction of RORγt and RORα (11, 41), #RORγt down of STAT3? induction by what?"
"5.7.3	Recent studies suggest the TGF-β promotes the expression of RORγt but represses its function (43). Only the additional presence of IL-6 or IL-21 signaling relieves the repression of RORγt and promotes the Th17 transcriptional program (29, 43)."
"29.6.1	In mice, Il23r is induced by IL-6 or IL-21 but is inhibited by high concentrations of TGF-β43. In human cells, IL23R expression was induced to some extent by IL-23 alone but not by IL-1b, consistent with a published report28 (Fig. 4e).",
"Simplify the addition of experimental data fitting the model: Transform an indirect relationship table to a direct one with intermediates, ensuring that the final model accout for all the interventions and their result"
]

# print input_xp_results([il_23,il_21,il_6,tgfb],[il_23r])
res=u"""
IL-23	IL-21	IL-6	TGF-β	IL-23R
0	0	0	0	0
0	1	0	0	0
0	0	0	1	0
0	0	0	2	0
0	0	0	3	0
0	0	1	3	0
0	1	0	3	0
0	0	1	0	1
0	0	1	2	1
0	1	0	1	1
0	1	0	2	1
0	0	1	1	2
1	1	0	1	2
1	1	0	2	2
1	0	1	1	2
1	0	1	2	2
1	0	1	3	1
1	1	0	3	0
"""

Experiment("21.4.5",[il_23,il_21,il_6,tgfb],[il_23r],res)

# print input_xp_results([il_23,il_21,il_6,tgfb],[il_17])
res="""
IL-23	IL-21	IL-6	TGF-β	IL-17
0	0	1	1	0
0	1	0	1	0
0	0	1	2	1
0	1	0	2	1
0	1	0	3	1
1	0	1	1	1
1	1	0	1	1
1	1	0	2	2
1	1	0	3	2
1	0	1	2	2
0	0	1	3	2
1	0	1	3	2
"""
Experiment("21.4.5",[il_23,il_21,il_6,tgfb],[il_17],res)

build_nx_graph(False)