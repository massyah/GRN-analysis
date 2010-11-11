#we need to generate all the possible state values of each input
""""
Simple by using itertools product 
http://docs.python.org/library/itertools.html#itertools.product

Do not require the graph, the influencing nodes are given in the focal tables
"""
import itertools

def domain(var):
	if var in ["IFNg","IFNgR","STAT1","Tbet"]:
		return 2
	return 1

def domain_for_vars(*args):
	ranges=[range(domain(x)+1) for x in args]
	vals=[x for x in itertools.product(*ranges)]
	return vals


def parse_focal_to_dict(ogene,matrix):
	lines=matrix.strip().split("\n")
	focal={}
	l=lines[0]
	assert(l[0].isalpha())
	#the two mappings
	idx_var=l.split("\t")[:-1]
	var_idx=dict(zip(idx_var,range(len(idx_var))))
	#generate all possible input values
	influencing=idx_var
	possible_input=domain_for_vars(*influencing)
	for v in possible_input:
		focal[tuple(zip(idx_var,v))]=0
	for l in lines[1:]:
		l=l.strip()
		if len(l)<1:
			continue
		row={}
		vals=l.split("\t")
		input_tuples=tuple(zip(idx_var,map(int,vals[:-1])))
		focal[tuple(input_tuples)]=int(vals[-1])
	return focal

focalParamsValues={}

def no_activator(tgt):
	return parse_focal_to_dict(tgt,"""
	K
	0
	""")

def one_activator(tgt,activator):
	return parse_focal_to_dict(tgt,"""
%s	K
1	1
"""%(activator))


focalParamsValues={}


focalParamsValues["IL18"]=no_activator("IL18")
focalParamsValues["IL12"]=no_activator("IL12")
focalParamsValues["IFNb"]=no_activator("IFNb")
focalParamsValues["IFNbR"]=one_activator("IFNbR","IFNb")
focalParamsValues["STAT6"]=one_activator("STAT6","IL4R")
focalParamsValues["IRAK"]=one_activator("IRAK","IL18R")


focalParamsValues["IL18R"]=parse_focal_to_dict("IL18R","""
IL18	STAT6	K
1	0	1
""")

focalParamsValues["IL4"]=parse_focal_to_dict("IL4","""
GATA3	STAT1	K
1	0	1
""")
focalParamsValues["IL4R"]=parse_focal_to_dict("IL4R","""
IL4	SOCS1	K
1	0	1
""")

focalParamsValues["STAT4"]=parse_focal_to_dict("STAT4","""
IL12R	GATA3	K
1	0	1
""")

focalParamsValues["GATA3"]=parse_focal_to_dict("GATA3","""
STAT6	Tbet	K
1	0	1
""")

focalParamsValues["IL12R"]=parse_focal_to_dict("IL12R","""
STAT6	IL12	K
0	1	1
""")

focalParamsValues["IFNg"]=parse_focal_to_dict("IFNg","""
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

focalParamsValues["IFNgR"]=parse_focal_to_dict("IFNgR","""
IFNg	SOCS1	K
1	0	1
1	1	1
2	1	1
2	0	2
""")

focalParamsValues["STAT1"]=parse_focal_to_dict("STAT1","""
IFNbR	IFNgR	K
1	0	1
0	1	1
1	1	1
0	2	2
1	2	2
""")


focalParamsValues["Tbet"]=parse_focal_to_dict("Tbet","""
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

focalParamsValues["SOCS1"]=parse_focal_to_dict("SOCS1","""
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


focalParamsValues["IFNg"]=parse_focal_to_dict("IFNg","""
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


def toMathematica():
	allFocals=[]
	for k,v in focalParamsValues.items():
		focals=[]
		for row,focal in v.items():
			inp=map(lambda x:"%s -> %i"%(x[0],x[1]),row)
			oup="%s -> %i"%(k,focal)
			focals.append("{"+",".join(inp)+"}->{"+oup+"}")
		allFocals.append("{"+",\n".join(focals)+"}")
	print "{"+",\n".join(allFocals)+"}"


ssTable=map(
	lambda gene:map(lambda k:k[0]+((gene,k[1]),),focalParamsValues[gene].items()),
	focalParamsValues.keys()
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
nodes=[[row for row in table if not hasTwoValue(row)] for table in ssTable]
def complement_state(ssState,depth):
	if depth==len(nodes): #finished
		ssState_map={}
		for k,v in ssState:
			if v!=0:
				ssState_map[k]=v
		print ssState_map
		return
	for newState in nodes[depth]:
		if ssAreCompatible(ssState,newState):
			complement_state(ssState+newState,depth+1)

complement_state(tuple(),0)