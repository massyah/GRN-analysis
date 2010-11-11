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
		self._values=[Any]*len(genes)
		if state:
				self._values=copy.copy(state._values)
	def __setitem__(self, key, value):
		self._values[key]=value
	def __getitem__(self,key):
		return self._values[key]
	def __mul__(self, other):
		temp=[Any]*N
		for i in geneRange:
			v1,v2=self._values[i],other._values[i]
			if v1!=v2:
				if v1!=Any and v2!=Any:
					return None
				elif v1==Any:
					temp[i]=v2
				elif v2==Any:
					temp[i]=v1
				else:
					assert False
			else:
				temp[i]=v1
		s=State()
		s._values=temp
		return s
	def __str__(self):
		return "s["+" ".join(map(str,self._values))+"]"
	def __hash__(self):
		return hash(tuple(self._values))
	def __eq__(self,other):
		return self._values==other._values
	
	
		
	
ft=reducedFocalParamsValues
# ft=bigFocalParamsValues

ssTable=map(lambda gene:map(lambda k:k[0]+((gene,k[1]),), ft[gene].items()), genes)



#we remove all ssStates where an item has two value mapped (the case for autoactivation)
ssNodes=[[row for row in table if not hasTwoValue(row)] for table in ssTable]
ssNodes.sort(key=len)

reducedNodes=[]
for table in ssNodes:
	rows=[]
	for state in table:
		s=State()
		for var,val in state:
			s[genesIndex[var]]=val
		rows.append(s)
	reducedNodes.append(rows)
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
	return results

# profilerun('compute_steady_state()')
results=compute_steady_state()
# print len(results),"steady states"
# for r in results:
# 	print r
