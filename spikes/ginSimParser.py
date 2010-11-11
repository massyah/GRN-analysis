import re,os,sys,itertools
from timeit import Timer
import networkx as nx
from IPython.Debugger import Tracer; debug_here = Tracer()

# inp="../2_nodes_model.ginml"
# inp="../ginSimModels/boolean_cell_cycle.ginml"
# inp="../ginSimModels/withBasal.ginml"
# inp="../ginSimModels/pairRule.ginml"
# inp="../ginSimModels/SP_1cell.ginml"
# inp="../ginSimModels/ErB2_model.ginml"
# inp="../ginSimModels/SP_6cells.ginml"
inp="../2006_mendoza.ginml"
# inp="../ginSimModels/TCRsig40.ginml"
# inp="../ginSimModels/Th_differentiation_full_annotated_model.ginml"


text=open(inp,"r").read()
#declare keys as variables
def generate_variables_for_keys(keyList):
	cmd=""
	for k in keyList:
		cmd+="%s=\"%s\"\n" %(k,k)
	return cmd

keyList="id basevalue maxvalue from_ to sign minvalue".split()
exec(compile(generate_variables_for_keys(keyList),"","exec"))

 
all_nodes_re=re.compile("""<node (.*?)>(.*?)</node>""",re.DOTALL)
all_edges_re=re.compile("""<edge (.*?)>(.*?)</edge>""",re.DOTALL)
interaction_param_re=re.compile("""<parameter (.*?)/>""",re.DOTALL)

def xml_param_to_tuples(m):
	key=m.group(1)
	if key in ["basevalue","minvalue","maxvalue","val"]: #ints
		val=int(m.group(2).strip("\""))
		return "(\"%s\",%i),"%(key,val)
	elif key in ["idActiveInteractions"]: #lists
		val=m.group(2).strip("\"").split(" ")
		return "(\"%s\",%s),"%(key,val.__repr__())
	else:
		val=m.group(2)
		return "(\"%s\",%s),"%(key,val)

nodes={}
edges={}

param_re=re.compile("(\w+)=(\"[^\"]+\")")

def xml_param_to_dict(param):
	tup=param_re.sub(xml_param_to_tuples,param)
	tup=dict(eval("("+tup+")"))
	return tup


#build a dict with all edge and nodes params
for m in all_nodes_re.finditer(text): 
	ndict=xml_param_to_dict(m.group(1))
	ndict["parameter"]=[]
	#search for interaction parameters
	for p in interaction_param_re.finditer(m.group(2)):
		interaction=xml_param_to_dict(p.group(1))
		ndict["parameter"].append(interaction)
	if "basevalue" not in ndict:
		ndict["basevalue"]=0
	nodes[ndict["id"]]=ndict
nxG=nx.DiGraph()
for m in all_edges_re.finditer(text):
	tup=xml_param_to_dict(m.group(1))
	nxG.add_edge(tup["from"],tup["to"])
	edges[tup["id"]]=tup
	

def domain(var):
	return nodes[var]["maxvalue"]

def domain_for_vars(*args):
	ranges=[range(domain(x)+1) for x in args]
	vals=[x for x in itertools.product(*ranges)]
	return vals
	
	
def influencing_genes(v):
	return list(set([x["from"] for x in edges.values() if x["to"]==v]))
	
#build the focal tables
def focal_table(var):
	influencing=influencing_genes(var)
	influencing.sort()
	table={}
	possible_input=domain_for_vars(*influencing)
	for v in possible_input:
		table[tuple(zip(influencing,v))]=0
	table[tuple(zip(influencing,[0]*len(influencing)))]=nodes[var]["basevalue"]
	for param in nodes[var]["parameter"]:
		activeVars=set()
		validDomain={}
		for interaction in param["idActiveInteractions"]:
			e=edges[interaction]
			if "maxvalue" in e:
				validDomain[e["from"]]=range(e["minvalue"],e["maxvalue"]+1)
			else:
				validDomain[e["from"]]=range(e["minvalue"],nodes[e["from"]]["maxvalue"]+1)
			activeVars.add(e["from"])
		for zeroVar in set(influencing)-activeVars:
			validDomain[zeroVar]=range(1) #should it be the base value?
		inputvals=[x for x in itertools.product(*validDomain.values())]
		for x in inputvals:
			inpVal=zip(tuple(validDomain.keys()),x)
			inpVal.sort()
			table[tuple(inpVal)]=param["val"]
	return table

focalParamsValues={}
genes=nx.dfs_preorder(nxG)

# genes=nx.dfs_preorder(nxG)
# genes=['IL12', 'IL18', 'IFNb', 'IFNbR', 'IL12R', 'IL18R', 'IRAK', 'STAT1', 'IL4', 'IL4R', 'STAT6', 'GATA3', 'STAT4', 'IFNg', 'IFNgR', 'SOCS1', 'Tbet']
# genes=['IFNb','IFNbR','IL18','IL18R','IRAK','IL12','IL12R','STAT4','STAT1', 'IL4', 'IL4R','SOCS1', 'STAT6','IFNg','IFNgR','GATA3','Tbet']
# genes=['Wg_external', 'Hh_external', 'Fz1', 'Dsh1', 'Nkd1', 'Slp1', 'En1', 'Ci1', 'Ciact1', 'Ptc1', 'Pka1', 'Cirep1', 'Hh1', 'Wg1']
# genes=['IL4_e', 'IL10RA', 'IL10RB', 'IL21_e', 'IL2_e', 'IFNG_e', 'IFNGR2', 'IL27_e', 'IL12_e', 'CGC', 'IFNB_e', 'IL6_e', 'IL15_e', 'TGFB_e', 'IL2RB', 'IL6RA', 'IL27RA', 'GP130', 'IL10_e', 'IL23_e', 'APC', 'IL15RA', 'IFNGR1', 'CD28', 'IFNBR', 'IL15R', 'IL27R', 'IL6R', 'TCR', 'IKB', 'NFAT', 'IFNGR', 'STAT1', 'IRF1', 'IL12RB1', 'IL10R', 'IL12R', 'IL21R', 'NFKB', 'STAT3', 'IL21', 'IL23', 'IL23R', 'STAT4', 'TGFBR', 'SMAD3', 'RORGT', 'FOXP3', 'IL2RA', 'IL2R', 'STAT5', 'IL4RA', 'IL4R', 'STAT6', 'IL12RB2', 'IL17', 'IL2', 'TGFB', 'proliferation', 'IL10', 'GATA3', 'TBET', 'RUNX3', 'IFNG', 'IL4']


#for the SP6 model, immediate answer, optimal, obtained with node order heuristic applied on Fz1
# genes=['Fz1', 'Dsh1', 'Nkd1', 'Slp1', 'En1', 'Ci1', 'Ciact1', 'Wg1', 'Cirep1', 'Hh1', 'Ptc1', 'Pka1', 'Fz2', 'Dsh2', 'Nkd2', 'Slp2', 'En2', 'Ci2', 'Ciact2', 'Wg2', 'Cirep2', 'Hh2', 'Ptc2', 'Pka2', 'Fz3', 'Dsh3', 'Nkd3', 'Slp3', 'En3', 'Ci3', 'Ciact3', 'Wg3', 'Cirep3', 'Hh3', 'Ptc3', 'Pka3', 'Fz4', 'Dsh4', 'Nkd4', 'Slp4', 'En4', 'Ci4', 'Ciact4', 'Wg4', 'Cirep4', 'Hh4', 'Ptc4', 'Pka4', 'Fz5', 'Dsh5', 'Nkd5', 'Slp5', 'En5', 'Ci5', 'Ciact5', 'Wg5', 'Cirep5', 'Hh5', 'Fz6', 'Dsh6', 'Nkd6', 'Ptc5', 'Pka5', 'Slp6', 'En6', 'Ci6', 'Ciact6', 'Ptc6', 'Pka6', 'Cirep6', 'Hh6', 'Wg6']



# genes=[x["id"] for x in sorted(nodes.values(),key=lambda x:(len(x["parameter"]),genes.index))]

# genes=nx.dfs_postorder(nxG)
# genes=nxG.nodes()
# genes=[x[0] for x in sorted(nxG.degree().items(),key=lambda x:x[1])]
# genes=[x[0] for x in sorted(nxG.in_degree().items(),key=lambda x:x[1])]
# genes=[x[0] for x in sorted(nx.in_degree_centrality(nxG).items(),key=lambda x:x[1])]


# genes.sort(key=lambda x:len(influencing_genes(x)))

print ",".join(genes)

for x in genes:
	focalParamsValues[x]=focal_table(x)
#copied from full_focal_table.py

ssTable=map(
	lambda gene:map(lambda k:k[0]+((gene,k[1]),),focalParamsValues[gene].items()),
	genes
	)
	
	
def ssAreCompatible(ss1,ss2): #assume the shortest is ss2
	d1=dict(ss1)
	for k,v in ss2: 
		if k in d1 and d1[k]!=v:
			return False
	return True
	
def hasTwoValue(l):
	tempDict={}
	for k,v in l:
		if k in tempDict and tempDict[k]!=v:
			return True
		tempDict[k]=v
	return False

#we remove all ssStates where an item has two value mapped (the case for autoactivation)
ssNodes=[[row for row in table if not hasTwoValue(row)] for table in ssTable]

# ssNodes.sort(key=len)
print ",".join([x[0][-1][0] for x in ssNodes])
# for s in ssNodes:
# 	s.sort(key=lambda x:x[-1],reverse=True)

maxDepth=0
depthVisit=[0]*(len(ssNodes)+1)
ssStates=[]
nComparisons=0
def complement_state(ssState,depth):
	global depthVisit,nComparisons
	depthVisit[depth]+=1
	ssStates
	if depth==len(ssNodes): #finished
		ssState_map={}
		for k,v in ssState:
			if v!=0:
				ssState_map[k]=v
		print ssState_map
		ssStates.append(ssState_map)
		return
	currState=dict(ssState)
	for newState in ssNodes[depth]:
		nComparisons+=1
		# if ssAreCompatible(ssState,newState): # we inline the call to ssAreCompatible
		goToNext=False
		for k,v in newState: 
			if k in currState and currState[k]!=v:
				goToNext=True
				break
		if goToNext:
			continue
		complement_state(ssState+newState,depth+1)

complement_state(tuple(),0)
realSS=[{'TCRphos': 1, 'TCRbind': 1, 'CD45': 1, 'Fyn': 1, 'IkB': 1, 'TCRlig': 1},
{'IkB': 1, 'TCRbind': 1, 'TCRlig': 1},
{'CD45': 1, 'IkB': 1, 'PAGCsk': 1},
{'IkB': 1, 'PAGCsk': 1},
{'IkB': 1, 'CD8': 1, 'TCRbind': 1, 'TCRlig': 1},
{'CD45': 1, 'IkB': 1, 'CD8': 1, 'PAGCsk': 1},
{'IkB': 1, 'CD8': 1, 'PAGCsk': 1}
]
# assert len(ssStates)==len(realSS)

print "Comp",nComparisons