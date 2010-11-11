#we describe a small part of activations in the th17 model and output that to ginSim
#we try to reproduce mendoza model

project="Th17 differentiation"
termFile=project+"/articles_th17.txt"
pubmed_files=[]


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

Activates(gata3,gata3)
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

"""
add the focal values manually of mendoza 2006 to generate a model comparable with ginsim original one
The edges are defined both by the predicates and by the focal tables
we store in all_edge_ids the set of all edges defined by the focal tables
an edge id is of the form STAT4_IFNg_0, which is an edge from STAT4 -> IFNg, with an activation value of [1,Max]
while Tbet_IFNg_1 has an activation value of [2,Max]
while Tbet_IFNg_0 has an activation value of [1,1]
how do they encode an edge with an activation value between 1 and 3, and 3 and 5 ?
the activation value of an edge is not the last element of the id, but a separate value in the edge definition
For example, 
	<node id="G0" basevalue="0" maxvalue="5">
	<node id="G1" basevalue="0" maxvalue="1">
	  <parameter idActiveInteractions="G0_G1_1" val="0"/>
	  <parameter idActiveInteractions="G0_G1_0" val="1"/>
	<edge id="G0_G1_0" from="G0" to="G1" minvalue="1" maxvalue="3" sign="positive">
	<edge id="G0_G1_1" from="G0" to="G1" minvalue="3" maxvalue="5" sign="positive">

This corresponds to this focal table for G1
G0	G1
[1,3]	1
[3,5]	0
"""

all_edges_ids
def parse_focal(txt):
	lines=txt.split("\n")
	focal=[]
	for l in lines:
		row={}
		if l[0].isalpha():
			vars=l.split("\t")
			var_idx=dict(zip(vars,range(len(vars))))
			continue
		vals=l.split("\t")
		for i in range(len(vals)):
			row[var_idx[i]]=vals[i]
		focal.append(row)
			
def parse_focal(gene,matrix):
	lines=matrix.split("\n")
	focal=[]
	for l in lines:
		l=l.strip()
		if len(l)<1:
			continue
		row={}
		if l[0].isalpha():
			#the two mappings
			idx_var=l.split("\t")
			var_idx=dict(zip(idx_var,range(len(idx_var))))
			continue
		vals=l.split("\t")
		edge_id=[]
		for i in range(len(vals)-1):
			if vals[i]!="0":
				edge_id.append("%s_%s_%i"%(idx_var[i],gene,int(vals[i])-1))
		edge_id=" ".join(edge_id)
		print """<parameter idActiveInteractions="%s" val="%i"/>"""% (edge_id,int(vals[-1]))

k_ifng=parse_focal("IFNg","""
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


exit(1)

build_nx_graph(onlyPositive=False)

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
edgeTemplate="""<edge id="%(id)s" from="%(src)s" to="%(tgt)s" minvalue="%(activationLevel)i" sign="%(sign)s">\n</edge>\n"""

def basal(id):
	return 0

def maxval(id):
	if id in [t.name for t in ifng,ifngr,stat1,tbet]:
		return 2
	return 1


def non_obvious_influences():
	"""output the nodes having multiple inputs
	For them, we have to define the focal values for a combination of input
	"""
	for n in nxG.nodes():
		influencing=nxG.predecessors(n)
		if len(influencing)>1:
			ifString=[]
			for i in influencing:
				nodeType=type(predicates_with_uid(nxG[i][n]["uid"])[0])
				if nodeType in [Downregulates,Inhibits,Inactivates,Sequestrates, BindsRepress]:
					ifString.append("-"+i)
				else:
					ifString.append("+"+i)
			print n,"influenced by",",".join(ifString)


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


nodes=nxG.nodes(data=True)
nodeList=" ".join([x[0] for x in nodes])
nodeDefs=""
edgeDefs=""

for n in nodes:
	nodeDefs+=nodeTemplate%{"id":n[0],"basal":basal(n[0]),"max":maxval(n[0])}
	influencing=nxG.predecessors(n[0])
	if len(influencing):
		for i in influencing:
			nodeType=type(predicates_with_uid(nxG[i][n[0]]["uid"])[0])
			if nodeType in [Downregulates,Inhibits,Inactivates,Sequestrates,BindsRepress]:
				focalVal=0
			else:
				focalVal=1
			nodeDefs+=nodeParams%{
			"edgesId":"_".join([i,n[0]]),
			"focal":focalVal
			}
	nodeDefs+="\n</node>\n"
	

for e in nxG.edges(data=True):
	p=predicates_with_uid(e[2]["uid"])[0]
	if type(p) in [Activates,Upregulates,Binds]:
		f=add_activation
	elif type(p) in [Downregulates,Inhibits,Inactivates,Sequestrates,BindsRepress]:
		f=add_inhibition
	else:
		f=add_activation
	edgeDefs+=f(e)


f=open("generatedModel.ginml","w")
f.write(template%{
	"nodeList":nodeList,
	"nodeDefs":nodeDefs,
	"edgesDefs":edgeDefs
	})
f.close()