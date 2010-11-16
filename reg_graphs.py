import omg_interface as omgI
import networkx as nx

def shortest_path(src,tgt):
	p=nx.shortest_path(nxG,src.name,tgt.name)
	justif=[]
	print "+->+".join(p)
	#assoc with predicates
	for i in range(len(p)-1):
		print nxG[p[i]][p[i+1]]
		pred=predicates_with_uid(nxG[p[i]][p[i+1]]["uid"])[0]
		justif.append(sentence_with_evidence_uid(pred.evidenceSentence)[0].string)
	print "\n".join(justif)


def all_paths(g,src,tgt):
	plist=[]
	path=[] #current path
	node_stack=[src]
	withData=True
	iterators=[g.edges_iter(src,withData)]

	while(len(node_stack))>0:
		node=node_stack[-1]
		iterator=iterators[-1]
		try:
			e=iterator.next()
			dest=e[1]
			if dest==tgt:
				path.append(e)
				plist.append(list(path))
				path.pop()
			elif dest not in node_stack:
				path.append(e)
				node_stack.append(dest)
				iterators.append(g.edges_iter(dest,withData))
		except StopIteration: #end of iteration
			node_stack.pop()
			iterators.pop()
			if len(path) > 0:
				path.pop()
	return plist
	
def draw_path(p):
	#add predecessors of every end nodes encountered
	g=nx.DiGraph()
	g.add_edges_from(p)
	for e in p:
		pred=predicates_with_uid(e[2]["uid"])[0]
		print type(pred)
		if e[0].startswith("AND"): #TODO: Fragile
			for t in nxG.predecessors(e[0]):
				g.add_edge(t,e[0],e[2])
	omgI.draw_nx_graph(g)

#prefix sort and cluster
#is this equiv to build the spanning tree rooted at the src and to trim what is downstream of the target ?
#no, because I can have a DAG

def all_paths_in_a_g(paths):
	bigG=nx.DiGraph()
	for p in paths:
		for e in p:
			bigG.add_edge(e[0],e[1],e[2])
	draw_path(bigG.edges(data=True))