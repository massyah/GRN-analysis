from cProfile import run as profilerun
import copy
# genes=['IL18', 'IL18R', 'IRAK', 'IFNg', 'IFNgR', 'STAT1', 'IL4', 'IL4R', 'STAT6', 'IL12R', 'STAT4', 'GATA3', 'Tbet', 'SOCS1', 'IL12', 'IFNbR', 'IFNb']

N=len(genes)
genesIndex=dict(zip(genes,range(len(genes))))
Any=-1
geneRange=range(N)

print "having",N,"genes"

class State(object):
	"""docstring for State"""
	def __init__(self, state=None):
		super(State, self).__init__()
		self._values={}
		if state:
				self._values=copy.copy(state._values)
	def __mul__(self, other):
		temp=copy.copy(self._values)
		for k,v2 in other._values.items():
			if k not in self._values:
				temp[k]=v2
			elif self._values[k]!=v2:
				return None
		s=State()
		s._values=temp
		return s
	def __str__(self):
		return "s{"+str(self._values)+"}"
	def __hash__(self):
		st=tuple(sorted(tuple(self._values.items()),key=lambda x:genesIndex[x[0]]))
		h=hash(st)
		return h
	def __eq__(self,other):
		return self._values==other._values

# ft=reducedFocalParamsValues
ft=bigFocalParamsValues

ssTable=map(lambda gene:map(lambda k:k[0]+((gene,k[1]),), ft[gene].items()), genes)



#we remove all ssStates where an item has two value mapped (the case for autoactivation)
ssNodes=[[row for row in table if not hasTwoValue(row)] for table in ssTable]


reducedNodes=[]
for table in ssNodes:
	rows=[]
	for state in table:
		s=State()
		for var,val in state:
			s._values[var]=val
		rows.append(s)
	reducedNodes.append(rows)
	
reducedNodes.sort(key=len)
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
			for s in col:
				res=r*s
				# print r,s,res
				if res:
					newResults.add(res)
		results=newResults
		# print "-"*12,len(results),"states"
	print "#",len(results),"steady states"
	return results

profilerun('compute_steady_state()')
# compute_steady_state()
# print len(results),"steady states"
# for r in results:
# 	print r
