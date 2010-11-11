import re,os,sys,itertools,copy
from cProfile import run as profilerun
from timeit import Timer
import networkx as nx
from IPython.Debugger import Tracer; debug_here = Tracer()

# from nodeOrderHeuristic import get_best_order

# inp="../2_nodes_model.ginml"
# inp="../ginSimModels/boolean_cell_cycle.ginml"
# inp="../ginSimModels/withBasal.ginml"
# inp="../ginSimModels/pairRule.ginml"
# inp="../ginSimModels/SP_1cell.ginml"
# inp="../ginSimModels/ErB2_model.ginml"

# inp="../2006_mendoza.ginml"
inp="../2006_mendoza_variant_without_tbet_ifng.ginml"

# inp="../ginSimModels/TCRsig40.ginml"
# inp="../ginSimModels/Th_differentiation_full_annotated_model.ginml"
# inp="../ginSimModels/SP_6cells.ginml"
# inp="../3_nodes_model_mutants.ginml"
		
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

#get the genes with the smallest focal table

# genes=get_best_order(nxG)
genes=nx.dfs_preorder(nxG)
#for the SP6 model, immediate answer, optimal, obtained with node order heuristic applied on Fz1
# genes=['Fz1', 'Dsh1', 'Nkd1', 'Slp1', 'En1', 'Ci1', 'Ciact1', 'Wg1', 'Cirep1', 'Hh1', 'Ptc1', 'Pka1', 'Fz2', 'Dsh2', 'Nkd2', 'Slp2', 'En2', 'Ci2', 'Ciact2', 'Wg2', 'Cirep2', 'Hh2', 'Ptc2', 'Pka2', 'Fz3', 'Dsh3', 'Nkd3', 'Slp3', 'En3', 'Ci3', 'Ciact3', 'Wg3', 'Cirep3', 'Hh3', 'Ptc3', 'Pka3', 'Fz4', 'Dsh4', 'Nkd4', 'Slp4', 'En4', 'Ci4', 'Ciact4', 'Wg4', 'Cirep4', 'Hh4', 'Ptc4', 'Pka4', 'Fz5', 'Dsh5', 'Nkd5', 'Slp5', 'En5', 'Ci5', 'Ciact5', 'Wg5', 'Cirep5', 'Hh5', 'Fz6', 'Dsh6', 'Nkd6', 'Ptc5', 'Pka5', 'Slp6', 'En6', 'Ci6', 'Ciact6', 'Ptc6', 'Pka6', 'Cirep6', 'Hh6', 'Wg6']


print ",".join(genes)


	
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


def print_state(s):
	print " ".join([x[0]+"="+str(x[1]) for x in s.items()]),

def merge_states(d1,d2):
	for k,v in d2.items():
		if k in d1 and d1[k]!=v:
			if d1[k]=="*":
				d1[k]=v
			elif v=="*":
				d2[k]=d1[k]
			else:
				return None
	res=d1.items()+d2.items()
	return res


def steady_state(ssNodes):
	possibleStates=map(dict,ssNodes[0])
	print possibleStates
	for i in range(1,len(genes)):
		# print i,":",len(possibleStates),"*",len(ssNodes[i])
		sys.stdout.flush()
		newStates=[]
		for p in possibleStates:
			pd=dict(p)
			for param in ssNodes[i]:
				nextState=merge_states(pd,dict(param))
				if nextState:
					newStates.append(nextState)
		possibleStates=set()
		for x in newStates:
			x.sort()
			possibleStates.add(tuple(x))
		# print possibleStates
	print len(possibleStates),"st states"
	return possibleStates
	


reducedFocalParamsValues={}
bigFocalParamsValues={}

for x in genes:
	reducedFocalParamsValues[x]=reduce_focal_table(focal_table(x))
	bigFocalParamsValues[x]=focal_table(x)

focalParamsValues=bigFocalParamsValues

ssTable=map(lambda gene:map(lambda k:k[0]+((gene,k[1]),), focalParamsValues[gene].items()), genes)


#we remove all ssStates where an item has two value mapped (the case for autoactivation)
ssNodes=[[row for row in table if not hasTwoValue(row)] for table in ssTable]
ssNodes.sort(key=len)

# profilerun('steady_state(ssNodes)')
# steady_state(ssNodes)