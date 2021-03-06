from cProfile import run as profilerun

def merge_states(d1,d2):
	global ncall
	ncall+=1
	d1=list(d1)
	for i in geneRange:
		if d1[i]!=d2[i]:
			if d1[i]=="_":
				d1[i]=d2[i]
			elif d2[i]=="_":
				continue
			else:
				return None
	return tuple(d1)


def fuse_table(t1,t2):
	newRows=set()
	for l1 in t1:
		for l2 in t2:
			r=merge_states(l1,l2)
			if r:
				newRows.add(r)
	return list(newRows)




def table_distance(t1,t2):
	if (len(t1)==0) or (len(t2)==0):
		return float('inf'),None
	if (len(t1)>300) or (len(t2)>300):
		return float('inf'),None
	ft=fuse_table(t1,t2)
	return len(ft),ft
		
		
def findmostsimilar(tables):
	minDist=float('inf')
	minPair=(None,None)
	minRes=None
	for i in range(len(tables)):
		for j in range(i+1,len(tables)):
			# d=table_distance(tables[i],tables[j])
			# assert dist[i,j]==d
			# assert dist[j,i]==d
			d,res=dist[i,j]
			if d<minDist:
				minPair=(i,j)
				minDist=d
				minRes=res
	# print minPair,maxProx,dist[minPair[0],minPair[1]]
	return minPair[0],minPair[1],minRes
	
def print_table(t):
	 for y in t:
		  print ",".join(map(str,y)),len([val for val in y if val!="_"])

def step(tables):
	global dist
	if len(tables)<=1:
		return True
	i1,i2,res=findmostsimilar(tables)
	if i1==None :
		print "finished",len(tables[0]),"st states"
		print_table(tables[0])
		return True
	t1=tables[i1]
	t2=tables[i2]
	nt=res
	tables[i1]=nt
	tables[i2]=[]
	#update the matrix distance:
	for i in range(N):
		d,res=table_distance(nt,tables[i])
		dist[(i,i1)]=d,res
		dist[(i1,i)]=d,res
		dist[(i2,i)]=float('inf'),None
		dist[(i,i2)]=float('inf'),None
	return False
	# print len(nt)

def print_distance_matrix():
	print "\t","\t".join(map(str,geneRange))
	for i in geneRange:
		print i,"\t"*(i+2),
		for j in range(i+1,N):
			print dist[(i,j)][0],"\t",
		print '\n' 

def sstate():
	global tables
	ncall=0
	i=0
	while(not step(tables)):
		# step(tables)
		print i
		sys.stdout.flush()
		# print len(tables),
		# print map(len,tables)
		i+=1




def rebuild_tables():
	global tables,dist,genes,N,geneRange,geneIndex
	ft=reducedFocalParamsValues

	genes=ft.keys()
	N=len(genes)
	geneRange=range(N)
	genesIndex=dict(zip(genes,range(len(genes))))
	# ft=bigFocalParamsValues
	ssTable=map(lambda gene:map(lambda k:k[0]+((gene,k[1]),), ft[gene].items()), genes)
	#we remove all ssStates where an item has two value mapped (the case for autoactivation)
	ssNodes=[[row for row in table if not hasTwoValue(row)] for table in ssTable]
	# ssNodes.sort(key=len)
	reducedNodes=[]
	for table in ssNodes:
		rows=[]
		for state in table:
			s=["_"]*N
			for var,val in state:
				s[genesIndex[var]]=val
			rows.append(tuple(s))
		reducedNodes.append(rows)

	tables=copy.copy(reducedNodes)

	dist={}
	for i in range(N):
		for j in range(i+1,N):
			d=table_distance(tables[i],tables[j])
			dist[(i,j)]=d
			dist[(j,i)]=d

rebuild_tables()
sstate()
# print dist
# profilerun("sstate()")