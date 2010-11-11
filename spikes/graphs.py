import networkx as nx
import omg_interface as omg
#example edges
# Binding(stat3,il21,True,8)
# Activation(il_6,stat3)
# Activation([smad3,nfat],foxp3)
# Activation(tgfb,smad3) #very old result
# SignalingAxis([il_2,il_2r,jak3,stat5])
# Binding(stat5,il17a)
# Binding(stat3,il17a)
# Binding(stat5,foxp3prom)
# Activation(tgfb,foxp3) #from PMID 14676299

speciesList="""ifng,ifngr,stat1,socs1,il-4,il-4r,stat6,gata3,t-bet,il-12,il-12r,stat4,il-18,il-18r,irak,ifnb,ifnbr,il21,stat3,tgfb,il-1b, il-12, il-4, il-6, il-5, il-18, il-2, il-21, il-23,il-17,il-23r il-2r smad3,nfat,foxp3,stat5,jak3,foxp3prom,runx1 il17a,il17f,rorgt"""

G=nx.DiGraph()
species=[s.strip() for s in speciesList.split(",")]
edges=[
("stat3","il21","12"),
("il-6","stat3",12),
("smad3","foxp3",1),
("nfat","foxp3",3),
("tgfb","smad3",4),
("il-2","il-2r",5),
("il-2r","jak3",6),
("jak3","stat5",1),
("stat5","il17a",2),
("stat3","il17a",3),
("tgfb","foxp3",4)
]
# print species


G.add_nodes_from(species)
for e in edges:
	G.add_edge(e[0],e[1],{"uid":12})
# G.add_edges_from(edges)

#iterate the edges to draw the G in omg
# omg.draw_nx_graph(G)

#shortest path bet il-2 and stat5
p=nx.shortest_path(G,"il-2","stat5")
print G.edges(p,data=True)
for i in range(len(p)-1):
	print p[i],p[i+1]
