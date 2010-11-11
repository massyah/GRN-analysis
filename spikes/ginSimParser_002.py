import re,os,sys,itertools,copy
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

def reduce_focal_table(g,focal):
	#we take the rows != 1
	active=[x for x in focal.items() if x[1]!=0]
	influence=influencing_genes(g)
	#we associate to each var in the table, the value it can take to activate the tgt
	newTable={}
	for gene in influence:
		possibleValues=set(range(1,domain(gene)+1)) #we do not account for 0 here
		for a in active:
			d=dict(a[0])
			if d[gene] in possibleValues:
				possibleValues.remove(d[gene])
		for v in possibleValues:
			row=dict(zip(influence,"*"*len(influence)))
			row[gene]=v
			newTable[tuple(row.items())]=0
	if len(newTable)==0:
		return focal
	for k,v in active:
		newTable[k]=v
	basal=tuple(zip(influence,[0]*len(influence)))
	newTable[basal]=nodes[g]["basevalue"]
	return newTable

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

genes=nx.dfs_preorder(nxG)
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
	res=dict(d1.items()+d2.items())
	res["History"]=[]
	if "History" in d1:
		res["History"]=copy.deepcopy(d1["History"])
	res["History"].append(d2)
	return res


def steady_state():
	possibleStates=map(dict,ssNodes[0])
	print possibleStates
	for i in range(1,len(genes)):
		print i,":",len(possibleStates),"*",len(ssNodes[i])
		sys.stdout.flush()
		newStates=[]
		for p in possibleStates:
			for param in ssNodes[i]:
				newDict=dict(param)
				nextState=merge_states(p,newDict)
				# if i==4:
				# 	print "\t",p,newDict,nextState
				if nextState:
					newStates.append(nextState)
		possibleStates=newStates
		# print possibleStates
	return possibleStates


bigFocalParamsValues={}
for x in genes:
	focalParamsValues[x]=reduce_focal_table(x,focal_table(x))
	bigFocalParamsValues[x]=focal_table(x)
	# focalParamsValues[x]=focal_table(x)

ssTable=map(lambda gene:map(lambda k:k[0]+((gene,k[1]),), focalParamsValues[gene].items()), genes)

#we remove all ssStates where an item has two value mapped (the case for autoactivation)
ssNodes=[[row for row in table if not hasTwoValue(row)] for table in ssTable]
ssStates=steady_state()




for s in ssStates:
	print [(x[0],x[1]) for x in s.items() if x[0]!="History" and x[1]!=0]
print len(ssStates)
sys.exit(0)
realSS=[{'cdh1': 1, 'p27': 1, 'Rb': 1}]


assert len(ssStates)==len(realSS)
for s in realSS:
	assert s in ssStates
for s in ssStates:
	assert s in realSS

