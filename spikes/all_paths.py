import networkx as nx

g=nx.DiGraph()
g.add_edges_from([
("A","B",{"name":"a"}),
("B","C"),
("A","C"),
("C","D"),
("B","A"),
("B","D"),
("C","A"),
("A","D")
])

src="A"
tgt="D"
plist=[]
path=[] #current path
node_stack=[src]
iterators=[g.edges_iter(src,data=True)]

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
			iterators.append(g.edges_iter(dest,data=True))
	except StopIteration: #end of iteration
		node_stack.pop()
		iterators.pop()
		if len(path) > 0:
			path.pop()
print plist