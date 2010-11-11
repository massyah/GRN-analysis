"""
No more global variables between modules, everything is passed along in the functions
"""
from cProfile import run as profilerun
import numpy,copy

# from ginSimParser_003 import *
from IPython.Debugger import Tracer; debug_here = Tracer()

Any=-1
rejected=0
dups=0
visited_states=set()
genes=None

def hasTwoValue(l):
	tempDict={}
	for k,v in l:
		if k in tempDict and tempDict[k]!=v:
			return True
		tempDict[k]=v
	return False

def merge_states(d1o,d2):
	global rejected
	d1=list(d1o)
	for i in geneRange1: #we do not look at the first pos, which is the uid tuple
		v1,v2=d1[i],d2[i]
		if v1!=v2:
			if v1==Any:
				d1[i]=v2
			elif v2==Any:
				continue
			else:
				# rejected+=1
				return None
	d1[0]=d1o[0]+d2[0]
	return tuple(d1)




def fuse_table(t1,t2,maxThreshold=None):
	global dups
	newRows=[]
	newRows_states=set()
	for l1 in t1:
		for l2 in t2:
			r=merge_states(l1,l2)
			# if r and r[1:] in newRows_states:
			# 	dups+=1
			if r and r[1:] not in newRows_states:
				newRows.append(r)
				newRows_states.add(r[1:])
				# newRows_states.add(r)
				# visited_states.add(r[1:])
			if maxThreshold and len(newRows_states)>=maxThreshold:
				return None
	return newRows


def table_distance(t1,t2,maxThreshold=None):
	if (len(t1)==0) or (len(t2)==0):
		return float('inf'),None
	if (len(t1)>300) or (len(t2)>300):
		return float('inf'),None
	ft=fuse_table(t1,t2,maxThreshold)
	if maxThreshold and ft==None:
		return float('inf'),None
	return len(ft),ft
		
		
def findmostsimilar(dist,tables):
	minDist=float('inf')
	minPair=(None,None)
	minRes=None
	for i in range(len(tables)):
		for j in range(i+1,len(tables)):
			d,res=dist[i*N+j]
			if d<minDist:
				minPair=(i,j)
				minDist=d
				minRes=res
	return minPair[0],minPair[1],minRes
	
def print_table(t):
	 for y in t:
		  print ",".join(map(str,y)),len([val for val in y if val!=Any])

def state_to_dict(ss):
	sstate=[]
	if type(ss[0])==tuple:
		ss=ss[1:]
	
	for i in geneRange:
		if ss[i]!=0:
			sstate.append((genes[i],ss[i]))
	sstate.sort()
	return dict(sstate)
	
def print_st_states(t):
	for ss in t:
		if type(ss[0])==tuple:
			ss=ss[1:]
		print state_to_dict(ss)
		
def step(tables,dist):
	if len(tables)<=1:
		return True
	i1,i2,res=findmostsimilar(dist,tables)
	if i1==None :
		return True
	t1=tables[i1]
	t2=tables[i2]
	nt=res
	tables[i1]=nt
	tables[i2]=[]
	#update the matrix distance:
	for i in range(N):
		d,res=table_distance(nt,tables[i])
		dist[i*N+i1]=d,res
		dist[i1*N+i]=d,res
		dist[i2*N+i]=float('inf'),None
		dist[i*N+i2]=float('inf'),None
	return False

def sstate(tables,dist):
	i=0
	while(not step(tables,dist)):
		i+=1
	return tables



def tableForMutations(focalTables,mutations):
	ft=copy.copy(focalTables)
	for gene,mutInterval in mutations.items():
		newTable={}
		for parameter,wtValue in ft[gene].items():
			closestMutVal=mutInterval[0]
			for mutVal in mutInterval:
				if abs(wtValue - mutVal) < abs(wtValue - closestMutVal):
					closestMutVal=mutVal
			newTable[parameter]=closestMutVal
		ft[gene]=newTable
	return ft
	
def rebuild_tables(focalTable,mutations):
	global genes,N,geneRange,geneRange1
	ft=tableForMutations(focalTable,mutations)
	genes=ft.keys()
	genes.sort()
	N=len(genes)
	geneRange=range(N)
	geneRange1=range(1,N+1)
	genesIndex=dict(zip(genes,range(len(genes))))

	ssTable=map(lambda gene:map(lambda k:k[0]+((gene,k[1]),), ft[gene].items()), genes)
	#we remove all ssStates where an item has two value mapped (the case for autoactivation)
	ssNodes=[[row for row in table if not hasTwoValue(row)] for table in ssTable]

	reducedNodes=[]
	uid=0
	for table in ssNodes:
		rows=[]
		for state in table:
			s=[Any]*N
			for var,val in state:
				s[genesIndex[var]]=val
			s.insert(0,(uid,))
			uid+=1
			rows.append(tuple(s))
		reducedNodes.append(rows)
	return reducedNodes

def trim_tables(tables,maxThreshold):
	changed=True
	nChanged=0
	tooBig=numpy.zeros((N,N))
	while changed:
		changed=False
		for i in range(N):
			if tables[i]==[]:continue
			for j in range(i+1,N):
				if tables[j]==[]:continue
				if tooBig[i,j]!=0:continue
				d,res=table_distance(tables[i],tables[j],maxThreshold)
				if res!=None:
					nChanged+=1
					tables[i]=res
					tables[j]=[]
					changed=True
					tooBig[i:]=0
					tooBig[j:]=0
					tooBig[:,i]=0
					tooBig[:,j]=0
				else:
					tooBig[i,j]=1
	return tables

def rebuild_distance_matrix(tables):
	dist=[(float('inf'),None)]*(N*N)
	for i in range(N):
		for j in range(i+1,N):
			d,res=table_distance(tables[i],tables[j])
			dist[i*N+j]=d,res
			dist[j*N+i]=d,res
	return dist
				
def count_non_null_items_in_distance_matrix(dist):
	c=0
	for i in geneRange:
		for j in range(i+1,N):
			if dist[i*N+j]!=(float('inf'),None):
				c+=1
	return c

def compute_sstate(focal_table, mutations):
	tables=rebuild_tables(focal_table,mutations)
	tables=trim_tables(tables,5)
	dist=rebuild_distance_matrix(tables)
	tables=sstate(tables,dist)
	return tables
	
# profilerun('compute_sstate()')
# compute_sstate()