#!/usr/bin/env python
# encoding: utf-8

"""
Model combining litterature mined results, ginsim output and stable states analysis

"""
import os,itertools,unicodedata,sys

project="Th17 differentiation"
termFile=project+"/articles_th17.txt"
pubmed_files=[f for f in os.listdir("./"+project) if f.startswith("entrez ")]
pubmed_files=[]

model_name="th17.ginml"

from IPython.Debugger import Tracer; debug_here = Tracer()
17

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
all_edges_id=[]
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

template="""<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE gxl SYSTEM "http://gin.univ-mrs.fr/GINsim/GINML_2_1.dtd">
<gxl xmlns:xlink="http://www.w3.org/1999/xlink">
	<graph id="default_name" class="regulatory" nodeorder="%(nodeList)s">
	%(nodeDefs)s
	%(edgesDefs)s
	</graph>
	</gxl>"""
nodeTemplate="""<node id="%(id)s" name="%(id)s" basevalue="%(basal)i" maxvalue="%(max)i">
"""
nodeParams="""<parameter idActiveInteractions="%(edgesId)s" val="%(focal)i"/>"""
edgeTemplate="""<edge id="%(id)s" from="%(src)s" to="%(tgt)s" minvalue="%(activationLevel)i" sign="%(sign)s">"""


def domain(var):
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
		print "\t".join(map(unicode,r))

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
	focalTxt=""
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
		edge_id=[]
		for i in range(len(vals)-1):
			row[idx_var[i]]=int(vals[i])
			if vals[i]!="0":
				edge_id.append(get_edge_id(idx_var[i],int(vals[i]),gene))
		if edge_id!=[]:
			edge_id=" ".join(edge_id)
			focalTxt+="""<parameter idActiveInteractions="%s" val="%i"/>\n"""% (edge_id,int(vals[-1]))
		else:
			focalTxt=""
				
		row=row.items()
		row.sort()
		focal[tuple(row)]=int(vals[-1])
	return (focal,focalTxt)




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



def add_activation(e):
	return edgeTemplate%{
	"id":"_".join(e[:2]),
	"src":e[0],
	"tgt":e[1],
	"activationLevel":1,
	"sign":"positive"
	}
def add_inhibition(e):
	return edgeTemplate%{
	"id":"_".join(e[:2]),
	"src":e[0],
	"tgt":e[1],
	"activationLevel":1,
	"sign":"negative"
	}


def get_edge_id(src,v1,tgt):
	global all_edges_id
	n_existing_edges=0
	for e in all_edges_id:
		if e[:3]==(src,v1,tgt):
			return e[3]
		elif (e[0],e[2])==(src,tgt):
			n_existing_edges+=1
	new_edge=(src,v1,tgt,"%s_%s_%i"%(src,tgt,n_existing_edges))
	all_edges_id.append(new_edge)
	return new_edge[3]

def build_ginSim_model():
	global all_edges_id
	nodeDefs=""
	edgeDefs=""
	nodes=nxG.nodes(data=True)
	nodeList=" ".join([x[0] for x in nodes])


	for n in nodes:
		gsName=ginSimNames(n[0])
		nodeDefs+=nodeTemplate%{"id":gsName,"basal":base_value(n[0]),"max":domain(n[0])}
		influencing=nxG.predecessors(n[0])
		if len(influencing):
			for i in influencing:
				nodeType=type(predicates_with_uid(nxG[i][n[0]]["uid"])[0])
				if n[0] in focalParams:
					nodeDefs+=focalParams[n[0]][1]
		if gsName in nodePositions:
			nodeDefs+="<nodevisualsetting>"+nodePositions[gsName]+"</nodevisualsetting>"
		nodeDefs+="\n</node>\n"
	
	edges={}
	for e in all_edges_id:
		etempl="""<edge id="%s" from="%s" to="%s" minvalue="%i" maxvalue="%i">"""%(e[3],ginSimNames(e[0]),ginSimNames(e[2]),e[1],e[1]) #TODO check if e[1]==max node value?
		if (e[0],e[2]) not in edges:
			edges[(e[0],e[2])]=[]
		edges[(e[0],e[2])].append(etempl)

	#add the inhibitions not appearing in the focal tables
	for e in nxG.edges(data=True):
		p=predicates_with_uid(e[2]["uid"])[0]
		if type(p) in [Downregulates,Inhibits,Inactivates,Sequestrates,BindsRepress]:
			if (e[0],e[1]) not in edges: #Weak: Assuming only one inhib between two nodes, real enough?
				edges[(e[0],e[1])]=[add_inhibition(e)]
			else:
				newDefs=[]
				for edg in edges[(e[0],e[1])]:
					newDefs.append(edg[:-1]+" sign=\"negative\">")
				edges[(e[0],e[1])]=newDefs
				print (e[0],e[1]),edges[(e[0],e[1])]
		else: #pos edge, add pos label to it
			newDefs=[]
			for edg in edges[(e[0],e[1])]:
				newDefs.append(edg[:-1]+" sign=\"positive\">")
			edges[(e[0],e[1])]=newDefs
			print (e[0],e[1]),edges[(e[0],e[1])]

	#add the edges visuals
	edgeDefs=""
	for e,v in edges.items():
		for anEdge in v:
			if e in edgePositions:
				edgeDefs+=anEdge+"\n<edgevisualsetting>\n"+edgePos[e]+"\n</edgevisualsetting>\n</edge>\n"
			else:
				edgeDefs+=anEdge+"</edge>\n"
				
	#build the full model
	model=template%{"nodeList":" ".join(map(ginSimNames,nxG.nodes())),
				"nodeDefs":nodeDefs,
				"edgesDefs":edgeDefs}
	#remove extended unicode chars to play well with gin sim
	model=unicodedata.normalize('NFKC', model).encode('ascii','ignore')
	# model=model.encode("utf-8")
	f=open(model_name,"w")
	f.write(model)
	f.close()

def verify_missing_focals():
	global nxG
	for gene in nxG.nodes():
		assert gene in focalParams

build_nx_graph(onlyPositive=False)

#visual defs, extracted with extractNodePos.py
def ginSimNames(n):
	replaceDict={tgfb.name:"TGF-b",rorgt.name:"ROR-gt"}
	if n in replaceDict:
		return replaceDict[n]
	else:
		return n

nodePositions={'FOXP3': '<point x="217" y="180"/>', 'Jak2-Jak1': '<point x="115" y="93"/>', 'STAT-3': '<point x="117" y="135"/>', 'IL-17': '<point x="152" y="238"/>', 'IL-6': '<point x="113" y="5"/>', 'ROR-gt': '<point x="117" y="184"/>', 'IL-21': '<point x="6" y="7"/>', 'IL-21R': '<point x="7" y="49"/>', 'Jak1-Jak3': '<point x="10" y="95"/>', 'IL-6R': '<point x="115" y="52"/>', 'TGF-b': '<point x="215" y="6"/>'}

edgePositions={}

#model definition & ginSim parameters
#from x-devonthink-item://207FB4BC-D2A1-4A52-A078-1086720A918F
focalParams={}


SignalingAxis([il_6,il_6r,jak2_jak1,stat3])
build_nx_graph(onlyPositive=False)
focalParams[il_6.name]=no_activator(il_6.name)
focalParams[il_6r.name]=one_activator(il_6r.name,il_6.name)
focalParams[jak2_jak1.name]=one_activator(jak2_jak1.name,il_6r.name)

#il 21 axis
axis=[il_21,il_21r,jak1_jak3,stat3]
SignalingAxis(axis)
Activates(stat3,il_21)
build_nx_graph(onlyPositive=False)

focalParams[il_21.name]=one_activator(il_21.name,stat3.name)
focalParams[il_21r.name]=one_activator(il_21r.name,il_21.name)
focalParams[jak1_jak3.name]=one_activator(jak1_jak3.name,il_21r.name)


#they both activate stat3
focalParams[stat3.name]=parse_focal(stat3.name,"""
Jak1-Jak3	Jak2-Jak1	K
0	0	0
0	1	1
1	0	1
1	1	1
""")




Activates(stat3, rorgt)
Inhibits(foxp3, rorgt)
build_nx_graph(onlyPositive=False)
focalParams[rorgt.name]=parse_focal(rorgt.name,"""
FOXP3	STAT-3	K
0	0	0
0	1	1
1	0	0
1	1	0
""")

Activates(tgfb, foxp3)
build_nx_graph(onlyPositive=False)
focalParams[foxp3.name]=one_activator(foxp3.name,tgfb.name)


Activates(rorgt, il_17)
build_nx_graph(onlyPositive=False)
focalParams[il_17.name]=one_activator(il_17.name,rorgt.name)

focalParams[tgfb.name]=no_activator(tgfb.name)



verify_missing_focals()

import fuse_tables as ft

#build the tables and compute the st states
focalTable,rft={},{}
for k,v in focalParams.items():
	focalTable[k]=focalParams[k][0]
genes=nxG.nodes()

for x in genes:
	rft[x]=reduce_focal_table(focalTable[x])

print rft
print ft.compute_sstate(rft,{})

#minimal req: 3 steady states, RORgt+,IL17+; Th0; Foxp3+
#for this, I need positive feedback loops