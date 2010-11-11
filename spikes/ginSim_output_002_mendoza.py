#we describe a small part of activations in the th17 model and output that to ginSim
#we try to reproduce mendoza model
import sys
project="Th17 differentiation"
termFile=project+"/articles_th17.txt"
pubmed_files=[]

from IPython.Debugger import Tracer; debug_here = Tracer()



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


#got from Mendoza using extractNodePos.py
nodePositions={'GATA3': '<ellipse x="266" y="317" width="75" height="30" backgroundColor="#FFCCCC" foregroundColor="#000000"/>',
 'IFNb': '<ellipse x="13" y="25" width="75" height="30" backgroundColor="#CC99FF" foregroundColor="#000000"/>',
 'IFNbR': '<ellipse x="14" y="133" width="75" height="30" backgroundColor="#99CCFF" foregroundColor="#000000"/>',
 'IFNg': '<rect x="140" y="25" width="75" height="30" backgroundColor="#CC99FF" foregroundColor="#000000"/>',
 'IFNgR': '<rect x="141" y="135" width="75" height="30" backgroundColor="#99CCFF" foregroundColor="#000000"/>',
 'IL12': '<ellipse x="401" y="27" width="75" height="30" backgroundColor="#CC99FF" foregroundColor="#000000"/>',
 'IL12R': '<ellipse x="402" y="138" width="75" height="30" backgroundColor="#99CCFF" foregroundColor="#000000"/>',
 'IL18': '<ellipse x="523" y="28" width="75" height="30" backgroundColor="#CC99FF" foregroundColor="#000000"/>',
 'IL18R': '<ellipse x="523" y="138" width="75" height="30" backgroundColor="#99CCFF" foregroundColor="#000000"/>',
 'IL4': '<ellipse x="261" y="25" width="75" height="30" backgroundColor="#CC99FF" foregroundColor="#000000"/>',
 'IL4R': '<ellipse x="264" y="136" width="75" height="30" backgroundColor="#99CCFF" foregroundColor="#000000"/>',
 'IRAK': '<ellipse x="523" y="228" width="75" height="30" backgroundColor="#99FF99" foregroundColor="#000000"/>',
 'SOCS1': '<ellipse x="13" y="226" width="75" height="30" backgroundColor="#99FF99" foregroundColor="#000000"/>',
 'STAT1': '<rect x="142" y="227" width="75" height="30" backgroundColor="#99FF99" foregroundColor="#000000"/>',
 'STAT4': '<ellipse x="404" y="229" width="75" height="30" backgroundColor="#99FF99" foregroundColor="#000000"/>',
 'STAT6': '<ellipse x="264" y="227" width="75" height="30" backgroundColor="#99FF99" foregroundColor="#000000"/>',
 'Tbet': '<rect x="146" y="316" width="75" height="30" backgroundColor="#FFCCCC" foregroundColor="#000000"/>'}

edgePos={('GATA3', 'IL4'): '<polyline points="303,332 365,328 366,41 298,40" line_style="straight" line_color="#009900" line_width="2" routage="manual"/>',
 ('GATA3', 'STAT4'): '<polyline points="303,332 445,339 441,244" line_style="straight" line_color="#FF0000" line_width="2" routage="auto"/>',
 ('GATA3', 'Tbet'): '<polyline points="303,332 245,326 183,331" line_style="straight" line_color="#FF0000" line_width="2" routage="manual"/>',
 ('IFNb', 'IFNbR'): '<polyline points="50,40 51,148" line_style="curve" line_color="#009900" line_width="2" routage="auto"/>',
 ('IFNbR', 'STAT1'): '<polyline points="51,148 116,146 116,227 179,242" line_style="straight" line_color="#009900" line_width="2" routage="auto"/>',
 ('IFNg', 'IFNgR'): '<polyline points="177,40 178,150" line_style="curve" line_color="#009900" line_width="2" routage="auto"/>',
 ('IFNgR', 'STAT1'): '<polyline points="178,150 179,242" line_style="curve" line_color="#009900" line_width="2" routage="manual"/>',
 ('IL12', 'IL12R'): '<polyline points="438,42 439,153" line_style="curve" line_color="#009900" line_width="2" routage="auto"/>',
 ('IL12R', 'STAT4'): '<polyline points="439,153 441,192 441,244" line_style="straight" line_color="#009900" line_width="2" routage="auto"/>',
 ('IL18', 'IL18R'): '<polyline points="560,43 560,153" line_style="curve" line_color="#009900" line_width="2" routage="auto"/>',
 ('IL18R', 'IRAK'): '<polyline points="560,153 566,197 560,243" line_style="straight" line_color="#009900" line_width="2" routage="auto"/>',
 ('IL4', 'IL4R'): '<polyline points="298,40 305,103 301,151" line_style="straight" line_color="#009900" line_width="2" routage="manual"/>',
 ('IL4R', 'STAT6'): '<polyline points="301,151 303,192 301,242" line_style="straight" line_color="#009900" line_width="2" routage="manual"/>',
 ('IRAK', 'IFNg'): '<polyline points="560,243 508,244 509,94 207,94 177,40" line_style="straight" line_color="#009900" line_width="2" routage="auto"/>',
 ('SOCS1', 'IFNgR'): '<polyline points="50,241 53,190 162,190 178,150" line_style="straight" line_color="#FF0000" line_width="2" routage="manual"/>',
 ('SOCS1', 'IL4R'): '<polyline points="50,241 68,205 283,205 301,151" line_style="straight" line_color="#FF0000" line_width="2" routage="manual"/>',
 ('STAT1', 'IL4'): '<polyline points="179,242 233,243 233,36 298,40" line_style="straight" line_color="#FF0000" line_width="2" routage="manual"/>',
 ('STAT1', 'SOCS1'): '<polyline points="179,242 122,243 50,241" line_style="straight" line_color="#009900" line_width="2" routage="manual"/>',
 ('STAT1', 'Tbet'): '<polyline points="179,242 167,292 183,331" line_style="straight" line_color="#009900" line_width="2" routage="manual"/>',
 ('STAT4', 'IFNg'): '<polyline points="441,244 401,244 400,111 188,111 177,40" line_style="straight" line_color="#009900" line_width="2" routage="auto"/>',
 ('STAT6', 'GATA3'): '<polyline points="301,242 286,299 303,332" line_style="straight" line_color="#009900" line_width="2" routage="manual"/>',
 ('STAT6', 'IL12R'): '<polyline points="301,242 319,187 424,186 439,153" line_style="straight" line_color="#FF0000" line_width="2" routage="auto"/>',
 ('STAT6', 'IL18R'): '<polyline points="301,242 333,204 548,204 560,153" line_style="straight" line_color="#FF0000" line_width="2" routage="manual"/>',
 ('Tbet', 'GATA3'): '<polyline points="183,331 241,342 303,332" line_style="straight" line_color="#FF0000" line_width="2" routage="manual"/>',
 ('Tbet', 'IFNg'): '<polyline points="183,331 104,328 102,41 177,40" line_style="straight" line_color="#009900" line_width="2" routage="auto"/>',
 ('Tbet', 'SOCS1'): '<polyline points="183,331 49,342 50,241" line_style="straight" line_color="#009900" line_width="2" routage="auto"/>',
 ('Tbet', 'Tbet'): '<polyline points="183,331 173,299 183,295 193,299 183,331" line_style="curve" line_color="#009900" line_width="2" routage="manual"/>'}



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

"""
add the focal values manually of mendoza 2006 to generate a model comparable with ginsim original one
The edges are defined both by the predicates and by the focal tables
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

we store in all_edge_ids the set of all edges defined by the focal tables, with the activation values, it is a mapping
	(src,v1,tgt)->edge_id
stored as the 4-uple src,v1,tgt,edge_id
	
"""

all_edges_id=[]
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
	focal=""
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
				edge_id.append(get_edge_id(idx_var[i],int(vals[i]),gene))
		edge_id=" ".join(edge_id)
		focal+="""<parameter idActiveInteractions="%s" val="%i"/>\n"""% (edge_id,int(vals[-1]))
	return focal
		


def one_activator(tgt,activator):
	return parse_focal(tgt,"""
%s	K
1	1
"""%(activator))


focalParams={}

focalParams["IL4"]=one_activator("IL4","GATA3")
focalParams["IL4R"]=one_activator("IL4R","IL4")
focalParams["IL12R"]=one_activator("IL12R","IL12")
focalParams["IL18R"]=one_activator("IL18R","IL18")
focalParams["STAT4"]=one_activator("STAT4","IL12R")
focalParams["IFNbR"]=one_activator("IFNbR","IFNb")
focalParams["STAT6"]=one_activator("STAT6","IL4R")
focalParams["IRAK"]=one_activator("IRAK","IL18R")
focalParams["GATA3"]=one_activator("GATA3","STAT6") #version without autoactivation



focalParams["IFNg"]=parse_focal("IFNg","""
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


focalParams["IFNgR"]=parse_focal("IFNgR","""
IFNg	SOCS1	K
1	0	1
1	1	1
2	1	1
2	0	2
""")

focalParams["STAT1"]=parse_focal("STAT1","""
IFNbR	IFNgR	K
1	0	1
0	1	1
1	1	1
0	2	2
1	2	2
""")


focalParams["Tbet"]=parse_focal("Tbet","""
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

focalParams["SOCS1"]=parse_focal("SOCS1","""
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

""" 
Unfortunately, in GinSim semantics, the sign of an edge is not uniquely defined, as it depends on the circuit that we consider. 
For this,we just consider that an interaction going to 0 is negative, but this is just not completely right.
//TODO: Check if we really need to output the max value
"""

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
edgeTemplate="""<edge id="%(id)s" from="%(src)s" to="%(tgt)s" minvalue="%(activationLevel)i" sign="%(sign)s">"""

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
			if n[0] in focalParams:
				nodeDefs+=focalParams[n[0]]
	if n[0] in nodePositions:
		nodeDefs+="<nodevisualsetting>"+nodePositions[n[0]]+"</nodevisualsetting>"
	nodeDefs+="\n</node>\n"
	

# for e in nxG.edges(data=True):
# 	p=predicates_with_uid(e[2]["uid"])[0]
# 	if type(p) in [Activates,Upregulates,Binds]:
# 		f=add_activation
# 	elif type(p) in [Downregulates,Inhibits,Inactivates,Sequestrates,BindsRepress]:
# 		f=add_inhibition
# 	else:
# 		f=add_activation
# 	edgeDefs+=f(e)


edges={}
for e in all_edges_id:
	etempl="""<edge id="%s" from="%s" to="%s" minvalue="%i" maxvalue="%i">"""%(e[3],e[0],e[2],e[1],e[1]) #TODO check if e[1]==max node value?
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
		if e in edgePos:
			edgeDefs+=anEdge+"\n<edgevisualsetting>\n"+edgePos[e]+"\n</edgevisualsetting>\n</edge>\n"
		else:
			edgeDefs+=anEdge+"</edge>\n"

f=open("1generatedModel.ginml","w")
f.write(template%{
	"nodeList":termList,
	"nodeDefs":nodeDefs,
	"edgesDefs":edgeDefs
	})
f.close()


