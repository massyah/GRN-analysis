"""We store a huge array of known incompatible states. To do so, each focal state is prefixed with a unique id, to perform a quick lookup in the incompat matrix
After merge, the merged state is incompatible with the union of the incompatible states. We can keep a list corresponding to the uids of each state where the merged state come from in the first position.
slower than without, and I only lower the number of calls to merge state by 3k, which is nothing.

In fact, I can only optimally save 25700 calls to merge states, ~ 20%, thereby lowering 2.6seconds to 2.08 secs

and there are only 1k dups

counting all unique visited states yields 85k states
"""
from cProfile import run as profilerun
import numpy,md5

Any=-1
rejected=0
dups=0
visited_states=set()
def merge_states(d1o,d2):
	global rejected
	d1=list(d1o)
	for i in range(1,N+1): #we do not look at the first pos, which is the uid tuple
		v1,v2=d1[i],d2[i]
		if v1!=v2:
			if v1==Any:
				d1[i]=v2
			elif v2==Any:
				continue
			else:
				rejected+=1
				return None
	d1[0]=d1o[0]+d2[0]
	return tuple(d1)




def fuse_table(t1,t2,maxThreshold=None):
	global dups
	newRows=[]
	newRows_states=set()
	for l1 in t1:
		for l2 in t2:
			# if len(l1[0])<5 and len(l2[0])<5 and incompatibleMatrix[l1[0],l2[0]].any() :
			# 	continue
			r=merge_states(l1,l2)
			if r and r[1:] in newRows_states:
				dups+=1
			if r and r[1:] not in newRows_states:
				newRows.append(r)
				newRows_states.add(r[1:])
				visited_states.add(r[1:])
			if maxThreshold and len(newRows)>=maxThreshold:
				return None
	return list(newRows)


def table_distance(t1,t2,maxThreshold=None):
	if (len(t1)==0) or (len(t2)==0):
		return float('inf'),None
	if (len(t1)>300) or (len(t2)>300):
		return float('inf'),None
	ft=fuse_table(t1,t2,maxThreshold)
	if maxThreshold and ft==None:
		return float('inf'),None
	return len(ft),ft
		
		
def findmostsimilar(tables):
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

def step(tables):
	global dist
	if len(tables)<=1:
		return True
	i1,i2,res=findmostsimilar(tables)
	if i1==None :
		print "finished",len(tables[0]),"st states"
		# print_table(tables[0])
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
	# print len(nt)

def print_distance_matrix():
	print "\t","\t".join(map(str,geneRange))
	for i in geneRange:
		print i,"\t"*(i+2),
		for j in range(i+1,N):
			print dist[i*N+j][0],"\t",
		print '\n' 

def sstate():
	global tables
	ncall=0
	i=0
	while(not step(tables)):
		# step(tables)
		print i,
		sys.stdout.flush()
		# print len(tables),
		# print map(len,tables)
		i+=1




def rebuild_tables():
	global tables,dist,genes,N,geneRange,geneIndex,incompatibleMatrix,allFocalStates
	ft=reducedFocalParamsValues
	genes=ft.keys()
	N=len(genes)
	geneRange=range(N)
	genesIndex=dict(zip(genes,range(len(genes))))

	ssTable=map(lambda gene:map(lambda k:k[0]+((gene,k[1]),), ft[gene].items()), genes)
	#we remove all ssStates where an item has two value mapped (the case for autoactivation)
	ssNodes=[[row for row in table if not hasTwoValue(row)] for table in ssTable]

	reducedNodes=[]
	uid=0
	totalFocalStates=sum(map(len,ssNodes))
	allFocalStates=[]
	geneToStates={}
	for table in ssNodes:
		rows=[]
		for state in table:
			s=[Any]*N
			for var,val in state:
				if var not in geneToStates:
					geneToStates[var]=set()
				geneToStates[var].add(uid)
				s[genesIndex[var]]=val
			s.insert(0,(uid,))
			uid+=1
			rows.append(tuple(s))
			allFocalStates.append(tuple(s))
		reducedNodes.append(rows)
	incompatibleMatrix=numpy.zeros((totalFocalStates,totalFocalStates))
	for k,v in geneToStates.items():
		v=list(v)
		for i in range(len(v)):
			for j in range(i+1,len(v)):
				s1,s2=v[i],v[j]
				if allFocalStates[s1][genesIndex[k]+1]!=allFocalStates[s2][genesIndex[k]+1]: #+1 du to the uid in the first pos
					incompatibleMatrix[s1,s2]=True
					incompatibleMatrix[s2,s1]=True
	tables=copy.copy(reducedNodes)
	print "tables built",len(tables)
	print map(len,tables)

def trim_tables(maxThreshold):
	global tables
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
	print "modified",nChanged
	print map(len,tables)

def rebuild_distance_matrix():
	global dist
	dist=[(float('inf'),None)]*(N*N)
	for i in range(N):
		for j in range(i+1,N):
			d,res=table_distance(tables[i],tables[j])
			dist[i*N+j]=d,res
			dist[j*N+i]=d,res

				
def count_non_null_items_in_distance_matrix():
	c=0
	for i in geneRange:
		for j in range(i+1,N):
			if dist[i*N+j]!=(float('inf'),None):
				c+=1
	return c

def compute_sstate():
	rebuild_tables()
	trim_tables(5)
	rebuild_distance_matrix()
	sstate()
profilerun('compute_sstate()')
# compute_sstate()