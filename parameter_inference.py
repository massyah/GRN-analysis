#!/usr/bin/env python
# encoding: utf-8

import os,itertools,unicodedata,sys

project="Th17 differentiation"
termFile=project+"/articles_th17.txt"
pubmed_files=[f for f in os.listdir("./"+project) if f.startswith("entrez ")]
pubmed_files=[]

model_name="paramInference_toy.ginml"

from IPython.Debugger import Tracer; debug_here = Tracer()


import omg_interface as omg
import reg_graphs as reg

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
#t IL-6


//Receptors
#t IL-6R


#t Jak2-Jak1
#t Jak1-Jak3

#t STAT3
"""

f=open("th17_model.terms","w")
f.write(termDefs.encode("utf-8"))
f.close()
parse_file("th17_model.terms")

cmd=generate_variable_command("""
jak2-jak1,jak1-jak3,stat3, il-6,il-6r
""")
exec(compile(cmd,"","exec"))

ranks={
"sign":[il_6],
"stats":[stat3],
"receptors":[il_6r]

}




experimentalResults=[]

def base_value(var):
	return 0

def domain(var):
	nonOne={}
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
			vals=map(int,l.split("\t"))
			res=tuple(zip(self.inpVars+self.outVars,vals))
			inpDict,outDict=dict(res[:len(self.inpVars)]),dict(res[len(self.inpVars):])
			self.relation.append((inpDict,outDict))
			all_results.append((inpDict,outDict,self))
	def influencing_nodes(self):
		global nxG
		res=[]
		res.extend(self.inpVars)
		for i in self.inpVars:
			for o in self.outVars:
				intermediates=reg.all_paths(nxG,i.name,o.name)
				for interPath in intermediates:
					for edge in interPath:
						res.append(allTerms[edge[0].lower()])
						res.append(allTerms[edge[1].lower()])
						#add the immediate predessors
						for pred in nxG.predecessors(edge[0]):
							res.append(allTerms[pred.lower()])
						for pred in nxG.predecessors(edge[1]):
							res.append(allTerms[pred.lower()])

		return list(set(res))

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

def parse_focal(gene,txt):
	global nxG,focalParams
	lines=[x for x in txt.split("\n") if len(x.strip())>0]
	focal={}

	if type(gene)==Term:
		gene=gene.name

	#fill with zero values
	input_vars=lines[0].strip().split("\t")[:-1]
	var_idx=dict(zip(input_vars,range(len(input_vars))))
	idx_var=lines[0].split("\t")

	possible_input=domain_for_vars(*input_vars)
	for v in possible_input:
		focal[tuple(zip(input_vars,v))]=0

	focal[tuple(zip(input_vars,[0]*len(input_vars)))]=base_value(gene)

	focalTxt=""
	for l in lines[1:]:
		l=l.strip()
		row={}
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
	focalParams[gene]=(focal,focalTxt)


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



		
focalParams={}

def no_activator(tgt):
	parse_focal(tgt,"""
K
0
	""")

def one_activator(tgt,activator):
	if type(activator)==Term:
		activator=activator.name
	parse_focal(tgt,"""
%s	K
1	1
"""%(activator))


#assumed interactions
Activates(il_6,jak2_jak1)
Activates(jak2_jak1,stat3)
Activates(jak1_jak3,stat3)





no_activator(il_6)
one_activator(jak2_jak1,il_6)
no_activator(jak1_jak3)

#we assume that stat3 is only activated when jak1_jak3 and jak2_jak1 are active
# parse_focal(stat3,u"""
# Jak1-Jak3	Jak2-Jak1	K
# 0	0	0
# 0	1	0 
# 1	0	0
# 1	1	1
# """)

#fitting

# print input_xp_results([il_6],[stat3])
res=u"""
IL-6	STAT3
0	0
1	0
1	1

"""
e1=Experiment("dummy",[il_6],[stat3],res)

#the problem is that we don't know the value of jak1_jak3 in the xp, we have to infer that jak1_jak3 must be active
build_nx_graph(False)


	
#we generate all the tuples for the influencing nodes that are in agreement with experimental result
influencing_nodes=e1.influencing_nodes()
influencing_nodes.sort(key=lambda x:x.name)
idx_to_name=dict(zip(range(len(influencing_nodes)),map(lambda x:x.name,influencing_nodes)))
name_to_idx=dict((v,k) for k, v in idx_to_name.iteritems())

rows_in_agreement=[]
for row in domain_for_vars(*influencing_nodes):
	for rel in e1.relation:
		valid=False
		for term in rel[0].keys(): #for the inputs
			if row[influencing_nodes.index(term)]!=rel[0][term]:
				valid=False
				break
			valid=True
		if not valid:
			continue
		for term in rel[1].keys(): #for the outputs
			if row[influencing_nodes.index(term)]!=rel[1][term]:
				valid=False
				break
			valid=True

		if valid: #we found a row in agreement with one row of the exp results
			rows_in_agreement.append(row)
print "rows in agreement"
for r in rows_in_agreement:
	print r
			
			
#we compare with the known focal tables
#we inverse each entry of the focal table
inv_focal={}
for k,v in focalParams.items():
	inv_map={}
	for k1, v1 in v[0].iteritems():
	    inv_map[v1] = inv_map.get(v1, [])
	    inv_map[v1].append(k1)
	inv_focal[k]=inv_map


#we build focal rows for each row in the influence_table
rows_in_agreement_with_focal=[]
for row in rows_in_agreement:
	row_focal={}
	for i in range(len(row)):
		var=idx_to_name[i]
		influencing_var_values=[]
		if len(nxG.predecessors(var))<1:
			continue
		for influencing_var in nxG.predecessors(var):
			influencing_var_values.append((influencing_var, row[name_to_idx[influencing_var]]))
		influencing_var_values.sort()
		row_focal[var]=[tuple(influencing_var_values),row[i]]
	#We compare the row_focal with the focalTable
	valid=True
	for var,focal_function in row_focal.items():
		if len(focal_function[0])==0: #it is a controllable input varialbe, we keep it
			continue
		if var not in focalParams:
			continue
		if focalParams[var][0][focal_function[0]]!=focal_function[1]:
			valid=False
			break
	if valid:
		print "keeping",row,"for ",row_focal
		rows_in_agreement_with_focal.append(row)
		
print "Rows in agreement with focal"
print rows_in_agreement_with_focal


import fuse_tables as ft

#build the tables and compute the st states
focalTable,rft={},{}
for k,v in focalParams.items():
	focalTable[k]=focalParams[k][0]
genes=nxG.nodes()

for x in genes:
	rft[x]=reduce_focal_table(focalTable[x])
	
tables=ft.compute_sstate(rft,{})
ft.print_st_states(tables[0])


#Prop checking
tables=ft.compute_sstate(rft,{il_6.name:[1,1]})
for state in tables[0]:
	print ft.state_to_dict(state)
	assert "STAT3" in ft.state_to_dict(state), "Test"
