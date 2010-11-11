from appscript import *
import networkx as nx
#by using omg edit->Copy as->Applescript menu item, it is very easy to build up some template code

omg=app("OmniGraffle Professional 5")
canv=app(u'OmniGraffle Professional 5').windows[1].canvas

def get_graphic(srcTerm):
	try:
		srcTerm.omg
	except:
		return srcTerm

	if srcTerm.omg==None:
		srcTerm.omg=canv.make(at=app.windows[1].canvas.graphics.end, new=k.shape, with_properties={k.text: {k.text: srcTerm.name, k.alignment: k.center}, k.draws_shadow: False, k.url: "http://www.google.com/search?q=%s"%(srcTerm.name),})
		srcTerm.omg.autosizing.set(k.full)
	return srcTerm.omg

def add_activation(srcTerm,tgtTerm):
	g1=get_graphic(srcTerm)
	g2=get_graphic(tgtTerm)
	lin=g1.connect(to=g2, with_properties={k.head_type: u'FilledArrow'})
	lin.url.set("http://www.google.com/search?q=%s %s"%(srcTerm.name,tgtTerm.name))
	
def add_upregulation(srcTerm,tgtTerm):
	g1=get_graphic(srcTerm)
	g2=get_graphic(tgtTerm)
	lin=g1.connect(to=g2, with_properties={k.head_type: u'FilledDoubleArrow'})
	lin.url.set("http://www.google.com/search?q=%s %s"%(srcTerm.name,tgtTerm.name))

	
def add_binding(srcTerm,tgtTerm):
	g1=get_graphic(srcTerm)
	g2=get_graphic(tgtTerm)
	lin=g1.connect(to=g2, with_properties={k.head_type: u'FilledCenterBall',k.head_scale:1.2})
	lin.url.set("http://www.google.com/search?q=%s %s"%(srcTerm.name,tgtTerm.name))
	
def add_and_node(predecessors):
	and_node=canv.make(at=app.windows[1].canvas.graphics.end, new=k.shape, with_properties={k.name:"Circle", k.draws_shadow: False,k.size:(14,14)})
	for src in predecessors:
		g=get_graphic(src)
		g.connect(to=and_node)
	return and_node

def add_inhibition(srcTerm,tgtTerm):
	g1=get_graphic(srcTerm)
	g2=get_graphic(tgtTerm)
	g1.connect(to=g2, with_properties={k.head_type: u'NegativeControls'})

def layout():
	canv.layout()
	
def draw_nx_graph(g):
	omg_mappings={}
	for e in g.edges_iter():
		print e
		if e[0] not in omg_mappings:
			omg_mappings[e[0]]=canv.make(at=app.windows[1].canvas.graphics.end, new=k.shape, with_properties={k.text: {k.text: e[0], k.alignment: k.center}, k.draws_shadow: False, k.url: "http://www.google.com/search?q=%s"%(e[0]),})
			omg_mappings[e[0]].autosizing.set(k.full)
		if e[1] not in omg_mappings:
			omg_mappings[e[1]]=canv.make(at=app.windows[1].canvas.graphics.end, new=k.shape, with_properties={k.text: {k.text: e[1], k.alignment: k.center}, k.draws_shadow: False, k.url: "http://www.google.com/search?q=%s"%(e[1]),})
			omg_mappings[e[1]].autosizing.set(k.full)
		omg_mappings[e[0]].connect(to=omg_mappings[e[1]], with_properties={k.head_type: u'FilledArrow'})
		if e[1].startswith("AND"):
			omg_mappings[e[1]].text.set("")
			omg_mappings[e[1]].name.set("Circle")
			omg_mappings[e[1]].size.set((14,14))
			omg_mappings[e[0]].autosizing.set(k.full)
	layout()