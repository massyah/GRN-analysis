from cProfile import run as profilerun
import copy
import nodeOrderHeuristic as nh

# genes=['IL18', 'IL18R', 'IRAK', 'IFNg', 'IFNgR', 'STAT1', 'IL4', 'IL4R', 'STAT6', 'IL12R', 'STAT4', 'GATA3', 'Tbet', 'SOCS1', 'IL12', 'IFNbR', 'IFNb']

#for the SP6 model, immediate answer, optimal, obtained with node order heuristic applied on Fz1
# genes=nh.get_best_order(nxG,"Fz1")
# genes=nh.get_best_order(nxG,"Slp1")
# genes=nh.get_best_order(nxG)
genes=nx.dfs_preorder(nxG)


N=len(genes)
genesIndex=dict(zip(genes,range(len(genes))))
Any="_"
geneRange=range(N)

print "having",N,"genes"

def merge_states(d1,d2):
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

ft=reducedFocalParamsValues
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

	
reducedNodes.sort(key=len)
# reducedNodes.sort(key=lambda x:sum(map(lambda row:len([val for val in row if val!="_"]),x)))

def print_table():
	for i in geneRange:
		print "table for",genes[i]
		for j in range(len(reducedNodes[i])):
			print i,j,reducedNodes[i][j]
		print "-"*12
#computing the steady states

def compute_steady_state():
	results=set(reducedNodes[0])
	for col in reducedNodes[1:]:
		newResults=set()
		for r in results:
			# print "   "+",".join(map(str,r))
			found=False
			for s in col:
				# print "* ",",".join(map(str,s)),
				res=merge_states(r,s)
				if res:
					# print "kept"
					newResults.add(res)
				else:
					# print 
					pass
		results=newResults
		# print "-"*12,len(results),"states"
	print "#",len(results),"steady states"
	return results

# profilerun('compute_steady_state()')

results=compute_steady_state()
print len(results),"steady states"
for r in results:
	print r
