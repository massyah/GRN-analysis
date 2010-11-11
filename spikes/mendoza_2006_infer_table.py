"""
We infer the focal table of Ifng that yields the desired phenotypes and mutant phenotypes
There 2*3*2*3=36 possible rows,
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

def blank_focal_table(gene):
	global nxG
	focal=[]
	#fill with zero values
	influencing=nxG.predecessors(gene)
	influencing.sort()
	possible_input=domain_for_vars(*influencing)
	focal.append(influencing+["K"])
	for v in possible_input:
		focal.append(v+(0,))
	# focal.append(tuple([0]*len(influencing)+[base_value(gene)]))
	focal.sort()
	for r in focal:
		print "\t".join(map(str,r))
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



# model without Tbet -> IFNg
# nxG.remove_edge("Tbet","IFNg")
# focalParamsValues["IFNg"]=parse_focal("IFNg","""
# IRAK	STAT4	K
# 0	0	0
# 0	1	0
# 1	0	0
# 1	1	1
# """)
# 

#model without STAT4 -> IFNg
# nxG.remove_edge("STAT4","IFNg")
# focalParamsValues["IFNg"]=parse_focal("IFNg","""
# IRAK	Tbet	K
# 0	1	1
# 0	2	2
# 1	1	1
# 1	2	2
# """)
#do not pass the IL-12 test

# #model without IFNb, IFNbr
# nxG.remove_edge("IFNbR","STAT1")
# focalParamsValues["STAT1"]=parse_focal("STAT1","""
# IFNgR	K
# 0	0
# 1	1
# 2	2
# """)
# #same as the full model

# #model without GATA3 -| STAT4
# nxG.remove_edge("GATA3","STAT4")
# focalParamsValues["STAT4"]=parse_focal("STAT4","""
# IL12R	K
# 0	0
# 1	1
# """)
# #same as the full model

#model without Tbet -| Gata3
# nxG.remove_edge("Tbet","GATA3")
# focalParamsValues["GATA3"]=parse_focal("GATA3","""
# STAT6	K
# 0	0
# 1	1
# """)
# #same as the full model
# nxG,focalParamsValues=remove_relation("Tbet","GATA3")

#model without STAT6 -| IL12R
# nxG.remove_edge("STAT6","IL12R")
# focalParamsValues["IL12R"]=parse_focal("IL12R","""
# IL12	K
# 0	0
# 1	1
# """)
# #fail the IL-12 mutant
# nxG,focalParamsValues=remove_relation("STAT6","IL12R")

def remove_relation(src,tgt):
	table=copy.deepcopy(focalParamsValues)
	graph=copy.deepcopy(nxG)
	graph.remove_edge(src,tgt)
	influencing=graph.predecessors(tgt)
	influencing.sort()
	possible_input=domain_for_vars(*influencing)
	fparams={}
	for v in possible_input:
		fparams[tuple(zip(influencing,v))]=0
	preserved=[]
	for interactions,value in table[tgt].items():
		inp=dict(interactions)
		if inp[src]!=0:
			continue
		del inp[src]
		newInteraction=inp.items()
		newInteraction.sort()
		fparams[tuple(newInteraction)]=value
	table[tgt]=fparams
	return graph,table



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
{'STAT6': 1, 'IL4': 1, 'GATA3': 1, 'IL4R': 1},
{},
{'IFNg': 2, 'IFNgR': 1, 'SOCS1': 1, 'Tbet': 2, 'STAT1': 1},
{'IFNg': 1, 'IFNgR': 1, 'SOCS1': 1, 'Tbet': 1, 'STAT1': 1}
])

il_12_p=({"IL12":[1,1]},[
{'STAT6': 1, 'IL12': 1, 'IL4R': 1, 'GATA3': 1, 'IL4': 1},
{'STAT4': 1, 'SOCS1': 1, 'STAT1': 1, 'IFNg': 2, 'IL12': 1, 'Tbet': 2, 'IFNgR': 1, 'IL12R': 1},
{'STAT4': 1, 'SOCS1': 1, 'STAT1': 1, 'IFNg': 1, 'IL12': 1, 'Tbet': 1, 'IFNgR': 1, 'IL12R': 1}
])

ifng_0=({"IFNg":[0,0]},[
{'STAT6': 1, 'IL4': 1, 'GATA3': 1, 'IL4R': 1},
{},
{'SOCS1': 1, 'Tbet': 2},
{'SOCS1': 1, 'Tbet': 1},
]
)


ifngr_0=({"IFNgR":[0,0]},[
{'STAT6': 1, 'IL4': 1, 'GATA3': 1, 'IL4R': 1},
{},
{'IFNg': 2, 'SOCS1': 1, 'Tbet': 2},
{'IFNg': 1, 'SOCS1': 1, 'Tbet': 1}

])

il_18_p=({"IL18":[1,1]},[
{'IL18': 1, 'STAT6': 1, 'IL4': 1, 'GATA3': 1, 'IL4R': 1},
{'IL18': 1, 'IRAK': 1, 'IL18R': 1},
{'IL18': 1, 'IRAK': 1, 'SOCS1': 1, 'STAT1': 1, 'IFNg': 2, 'Tbet': 2, 'IFNgR': 1, 'IL18R': 1},
{'IL18': 1, 'IRAK': 1, 'SOCS1': 1, 'STAT1': 1, 'IFNg': 1, 'Tbet': 1, 'IFNgR': 1, 'IL18R': 1}
])

il_12_18_p=({"IL12":[1,1],"IL18":[1,1]},[
{'IL18': 1, 'STAT6': 1, 'IL4R': 1, 'IL12': 1, 'IL4': 1, 'GATA3': 1},
{'IL18': 1, 'IRAK': 1, 'STAT4': 1, 'SOCS1': 1, 'STAT1': 1, 'IFNg': 2, 'IL12': 1, 'Tbet': 2, 'IFNgR': 1, 'IL18R': 1, 'IL12R': 1},
{'IL18': 1, 'IRAK': 1, 'STAT4': 1, 'SOCS1': 1, 'STAT1': 1, 'IFNg': 2, 'IL12': 1, 'Tbet': 1, 'IFNgR': 1, 'IL18R': 1, 'IL12R': 1}

])
gata3_p=({"GATA3":[1,1]},[
{'STAT6': 1, 'IL4': 1, 'GATA3': 1, 'IL4R': 1},
{'SOCS1': 1, 'STAT1': 1, 'IFNg': 2, 'Tbet': 2, 'IFNgR': 1, 'GATA3': 1},
{'SOCS1': 1, 'STAT1': 1, 'IFNg': 1, 'Tbet': 1, 'IFNgR': 1, 'GATA3': 1},
])


mutationsPhenotypes=[wt,il_12_p,ifng_0,ifngr_0,il_18_p,il_12_18_p,gata3_p]
# mutationsPhenotypes=[wt]
# relationShip_to_test=[("STAT6","IL12R"),("Tbet","GATA3")]
relationShip_to_test=nxG.edges()
for r in relationShip_to_test:
	nxGtemp,focalTemp=remove_relation(*r)
	valid=True

	for m in mutationsPhenotypes:
		tables=ft.compute_sstate(focalTemp,m[0])
		ssStates=map(ft.state_to_dict,tables[0])
	
		for ssState in ssStates:
			valid &= (ssState in m[1])


		for ssState in m[1]:
			valid&=(ssState in ssStates)
		if not valid:
			break

	if valid:
		print "relation",r,"not needed"
	else:
		print "relation",r,"needed"