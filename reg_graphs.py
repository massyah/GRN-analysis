from appscript import *
omg=app("OmniGraffle Professional 5")
canv=app(u'OmniGraffle Professional 5').windows[1].canvas


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
	draw_nx_graph(g)
	
def draw_data_less_nx_graph(g):
	omg_mappings={}
	for e in g.edges_iter(data=True):
		print e

		if e[0] not in omg_mappings:
			omg_mappings[e[0]]=canv.make(at=app.windows[1].canvas.graphics.end, new=k.shape, with_properties={k.text: {k.text: e[0], k.alignment: k.center}, k.draws_shadow: False, k.url: "http://www.google.com/search?q=%s"%(e[0]),})
			omg_mappings[e[0]].autosizing.set(k.full)
		if e[1] not in omg_mappings:
			omg_mappings[e[1]]=canv.make(at=app.windows[1].canvas.graphics.end, new=k.shape, with_properties={k.text: {k.text: e[1], k.alignment: k.center}, k.draws_shadow: False, k.url: "http://www.google.com/search?q=%s"%(e[1]),})
			omg_mappings[e[1]].autosizing.set(k.full)
		head_type=u'FilledArrow'
		omg_mappings[e[0]].connect(to=omg_mappings[e[1]], with_properties={k.head_type: head_type})

		if e[1].startswith("AND"):
			omg_mappings[e[1]].text.set("")
			omg_mappings[e[1]].name.set("Circle")
			omg_mappings[e[1]].size.set((14,14))
			omg_mappings[e[0]].autosizing.set(k.full)
	canv.layout()
def draw_nx_graph(g):
	omg_mappings={}
	for e in g.edges_iter(data=True):
		print e
		pred=predicates_with_uid(e[2]["uid"])[0]
		print e,type(pred)
		if e[0] not in omg_mappings:
			omg_mappings[e[0]]=canv.make(at=app.windows[1].canvas.graphics.end, new=k.shape, with_properties={k.text: {k.text: e[0], k.alignment: k.center}, k.draws_shadow: False, k.url: "http://www.google.com/search?q=%s"%(e[0]),})
			omg_mappings[e[0]].autosizing.set(k.full)
		if e[1] not in omg_mappings:
			omg_mappings[e[1]]=canv.make(at=app.windows[1].canvas.graphics.end, new=k.shape, with_properties={k.text: {k.text: e[1], k.alignment: k.center}, k.draws_shadow: False, k.url: "http://www.google.com/search?q=%s"%(e[1]),})
			omg_mappings[e[1]].autosizing.set(k.full)
		#relation type
		if type(pred) in [Activates, Binds, Complexes, Phosphorylates, RequiredToActivate, Upregulates]:
			head_type=u'FilledArrow'
		else:
			head_type=u'NegativeControls'
		omg_mappings[e[0]].connect(to=omg_mappings[e[1]], with_properties={k.head_type: head_type})

		if e[1].startswith("AND"):
			omg_mappings[e[1]].text.set("")
			omg_mappings[e[1]].name.set("Circle")
			omg_mappings[e[1]].size.set((14,14))
			omg_mappings[e[0]].autosizing.set(k.full)
	canv.layout()
	
#prefix sort and cluster
#is this equiv to build the spanning tree rooted at the src and to trim what is downstream of the target ?
#no, because I can have a DAG

def all_paths_in_a_g(paths):
	bigG=nx.DiGraph()
	for p in paths:
		for e in p:
			bigG.add_edge(e[0],e[1],e[2])
	draw_path(bigG.edges(data=True))