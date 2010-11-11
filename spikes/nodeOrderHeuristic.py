import random

def add_next(nxG,visited):
	possibleNext=[]
	for x in visited:
		for y in nxG.successors(x):
			possibleNext.append(y)
	# for x in possibleNext:
	requirements=[(pn,nxG.predecessors(pn)) for pn in set(possibleNext)-set(visited)]
	if len(requirements)==0:
		return None
	nextState=min([(len(set(r[1])-set(visited)),r[0]) for r in requirements])[1]
	return nextState


def get_best_order(nxG,start=None):

	tovisit=[x[0] for x in sorted(nxG.in_degree().items(),key=lambda x:x[1])]
	N=len(tovisit)
	visited=[]
	

	minDegree=min(nxG.in_degree().items(),key=lambda x:x[1])[1]
	
	
	#we take all the nodes with minimal degree
	if minDegree==0:
		while nxG.in_degree(tovisit[0])==0:
			visited.append(tovisit.pop(0))
	elif start:
		visited=[start]
		tovisit.remove(start)
	else:
		a=random.choice(tovisit)
		print "taking any",a
		visited.append(a)
		tovisit.remove(a)
		
	visitedp=None
	while len(visited)!=N:
		assert len(visited)>0
		# visited.append(tovisit.pop(0))
		while visitedp!=visited:
			visitedp=[x for x in visited] #fast copy
			nextGene=add_next(nxG,visited)
			if nextGene==None:
				break
			visited.append(nextGene)
			tovisit.remove(nextGene)

	assert len(visited)==len(nxG.nodes())
	assert len(tovisit)==0
	return visited
