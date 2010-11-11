"""
We automate the verification of the WT and mutant phenotypes by using the stables states
We define the model
"""
from cProfile import run as profilerun
import numpy
import sys
import re,os,sys,itertools,copy



project="Th17 differentiation"
termFile=project+"/articles_th17.txt"
pubmed_files=[]


#miner_007 should be run once before running this

from IPython.Debugger import Tracer; debug_here = Tracer()



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



termList="IFNg IFNgR STAT1 Tbet SOCS1 IFNb IFNbR IL18 IL18R IRAK IL12 IL12R STAT4 IL4 IL4R STAT6 GATA3"
terms=[Term(x,[]) for x in termList.split(" ")]


cmd=generate_variable_command(termList.lower(),sep=" ")
exec(compile(cmd,"","exec"))

#signaling axis from the 2006 paper
SignalingAxis([ifng,ifngr,stat1,socs1])
SignalingAxis([il4,il4r,stat6,gata3,il4])
SignalingAxis([il12,il12r,stat4,ifng])
SignalingAxis([il18,il18r,irak,ifng])
SignalingAxis([ifnb,ifnbr,stat1])

Activates(stat1,tbet)
Activates(tbet,ifng)

# Activates(gata3,gata3)
Activates(tbet,tbet)
Activates(tbet,socs1)



Inhibits(gata3,tbet)
Inhibits(tbet,gata3)
Inhibits(gata3,stat4)
Inhibits(socs1,ifngr)

Inhibits(stat1,il4)
Inhibits(socs1,il4r)
Inhibits(stat6,il12r)
Inhibits(stat6,il18r)

build_nx_graph(onlyPositive=False)

def domain(var):
	if var in ["IFNg","IFNgR","STAT1","Tbet"]:
		return 2
	else:
		return 1
		
def base_value(var):
	return 0
def domain_for_vars(*args):
	ranges=[range(domain(x)+1) for x in args]
	vals=[x for x in itertools.product(*ranges)]
	return vals


def reduce_focal_table(ft):
	genes=set()
	for row in ft.items():
		for col in row[0]:
			genes.add(col[0])

	inputs=set([x[0] for x in ft.items()])
	newRows={}
	for g in genes:
		for i in range(domain(g)+1):
			val=(g,i)
			focals=set([x[1] for x in ft.items() if val in x[0]])
			if len(focals)==1:
				v=focals.pop()
				inpToRemove=[x[0] for x in ft.items() if val in x[0]]
				inputs=inputs-set(inpToRemove)
				newRows[(val,)]=v
	for i in inputs:
		newRows[i]=ft[i]
	return newRows

def parse_focal(gene,txt):
	global nxG
	lines=txt.split("\n")
	focal={}
	#fill with zero values
	influencing=nxG.predecessors(gene)
	influencing.sort()
	possible_input=domain_for_vars(*influencing)
	for v in possible_input:
		focal[tuple(zip(influencing,v))]=0
	focal[tuple(zip(influencing,[0]*len(influencing)))]=base_value(gene)
	
	for l in lines:
		l=l.strip()
		if len(l)<1:
			continue
		
		row={}
		if l[0].isalpha():
			vars=l.split("\t")
			var_idx=dict(zip(vars,range(len(vars))))
			idx_var=l.split("\t")
			continue
		vals=l.split("\t")
		for i in range(len(vals)-1):
			row[idx_var[i]]=int(vals[i])
		row=row.items()
		row.sort()
		focal[tuple(row)]=int(vals[-1])
	return focal


focalParamsValues={}

def no_activator(tgt):
	return parse_focal(tgt,"""
	K
	0
	""")

def one_activator(tgt,activator):
	return parse_focal(tgt,"""
%s	K
1	1
"""%(activator))

genes=['IL18', 'IL18R', 'IRAK', 'IFNg', 'IFNgR', 'STAT1', 'IL4', 'IL4R', 'STAT6', 'IL12R', 'STAT4', 'GATA3', 'Tbet', 'SOCS1', 'IL12', 'IFNbR', 'IFNb']

focalParamsValues={}


focalParamsValues["IL18"]=no_activator("IL18")
focalParamsValues["IL12"]=no_activator("IL12")
focalParamsValues["IFNb"]=no_activator("IFNb")
focalParamsValues["IFNbR"]=one_activator("IFNbR","IFNb")
focalParamsValues["STAT6"]=one_activator("STAT6","IL4R")
focalParamsValues["IRAK"]=one_activator("IRAK","IL18R")


focalParamsValues["IL18R"]=parse_focal("IL18R","""
IL18	STAT6	K
1	0	1
""")

focalParamsValues["IL4"]=parse_focal("IL4","""
GATA3	STAT1	K
1	0	1
""")
focalParamsValues["IL4R"]=parse_focal("IL4R","""
IL4	SOCS1	K
1	0	1
""")

focalParamsValues["STAT4"]=parse_focal("STAT4","""
IL12R	GATA3	K
1	0	1
""")

focalParamsValues["GATA3"]=parse_focal("GATA3","""
STAT6	Tbet	K
1	0	1
""")

focalParamsValues["IL12R"]=parse_focal("IL12R","""
STAT6	IL12	K
0	1	1
""")

focalParamsValues["IFNg"]=parse_focal("IFNg","""
STAT4	Tbet	IRAK	K
1	0	0	1
0	1	0	1
1	1	0	1
0	2	0	2
0	2	1	2
1	2	1	2
1	0	1	2
1	1	1	2
1	1	1	2
0	1	1	1
1	2	0	2
""")

focalParamsValues["IFNgR"]=parse_focal("IFNgR","""
IFNg	SOCS1	K
1	0	1
1	1	1
2	1	1
2	0	2
""")

focalParamsValues["STAT1"]=parse_focal("STAT1","""
IFNbR	IFNgR	K
1	0	1
0	1	1
1	1	1
0	2	2
1	2	2
""")


focalParamsValues["Tbet"]=parse_focal("Tbet","""
STAT1	GATA3	Tbet	K
1	0	0	1
1	1	1	1
1	0	1	1
0	0	1	1
1	1	2	2
1	0	2	2
2	0	0	2
2	0	1	2
2	0	2	2
0	0	2	2
""")

focalParamsValues["SOCS1"]=parse_focal("SOCS1","""
Tbet	STAT1	K
1	0	1
0	1	1
1	1	1
1	2	1
2	1	1
2	0	1
0	2	1
2	2	1
""")


focalParamsValues["IFNg"]=parse_focal("IFNg","""
STAT4	Tbet	IRAK	K
1	0	0	1
0	1	0	1
1	1	0	1
0	2	0	2
0	2	1	2
1	2	1	2
1	0	1	2
1	1	1	2
1	1	1	2
0	1	1	1
1	2	0	2
""")

focalParamsValues["IFNg"]=parse_focal("IFNg","""
STAT4	Tbet	IRAK	K
1	0	0	1
0	1	0	1
1	1	0	1
0	2	0	2
0	2	1	2
1	2	1	2
1	0	1	2
1	1	1	2
1	1	1	2
0	1	1	1
1	2	0	2
""")


reducedFocalParamsValues={}
for x in genes:
	reducedFocalParamsValues[x]=reduce_focal_table(focalParamsValues[x])
	
N=len(genes)
genesIndex=dict(zip(genes,range(len(genes))))
Any=-1
geneRange=range(N)
geneRange1=range(N+1)

import fuse_tables_008 as ft


wt=({},[
(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
(0, 0, 2, 1, 0, 0, 0, 0, 1, 1, 2, 0, 0, 0, 0, 0, 0),
(0, 0, 1, 1, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0),
(0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0)
])
il_12_p=({"IL12":[1,1]},[
(0, 0, 2, 1, 0, 0, 0, 1, 1, 1, 2, 1, 0, 0, 1, 0, 0),
(0, 0, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 1, 0, 0),
(0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0)
])

ifng_0=({"IFNg":[0,0]},[
(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
(0, 0, 2, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0),
(0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0),
(0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0)
]
)


ifngr_0=({"IFNgR":[0,0]},[
(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
(0, 0, 2, 0, 0, 0, 0, 0, 1, 0, 2, 0, 0, 0, 0, 0, 0),
(0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0),
(0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0)
])

il_18_p=({"IL18":[1,1]},[
(0, 1, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0),
(1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1),
(1, 0, 2, 1, 0, 1, 0, 0, 1, 1, 2, 0, 0, 0, 0, 0, 1),
(1, 0, 1, 1, 0, 1, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 1)
])

il_12_18_p=({"IL12":[1,1],"IL18":[1,1]},[
(1, 0, 2, 1, 0, 1, 0, 1, 1, 1, 2, 1, 0, 0, 1, 0, 1),
(1, 0, 1, 1, 0, 1, 0, 1, 1, 1, 2, 1, 0, 0, 1, 0, 1),
(0, 1, 0, 0, 1, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0)
])
gata3_p=({"GATA3":[1,1]},[
(0, 0, 2, 1, 1, 0, 0, 0, 1, 1, 2, 0, 0, 0, 0, 0, 0),
(0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0),
(0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0)
])


mutationsPhenotypes=[wt,il_12_p,ifng_0,ifngr_0,il_18_p,il_12_18_p,gata3_p]

for m in mutationsPhenotypes:
	tables=ft.compute_sstate(reducedFocalParamsValues,m[0])
	print m[0],
	valid=True
	for ssState in tables[0]:
		valid = valid & (ssState[1:] in m[1])
	if not valid:
		print "diff on"
		ft.print_st_states(tables[0])
		print "not equiv to"
		ft.print_st_states(m[1])
	print valid,
	print "-"*8