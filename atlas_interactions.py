import re,os
import networkx as nx

atlas_file="/Volumes/Gauss/Documents/VDR - STAT3 interactions/mmc2_human.csv"
data=map(lambda x:x.split(";"),open(atlas_file).readlines()[3:])
#col 2 and 4 for the prot name


ppi_nxG=nx.Graph()
for l in data:
	ppi_nxG.add_edge(l[2],l[4])
	
print nx.shortest_path(ppi_nxG,"VDR","STAT3")
stat3_neigh=ppi_nxG.neighbors("STAT3")
vdr_neigh=ppi_nxG.neighbors("VDR")

print set(stat3_neigh).intersection(set(vdr_neigh)) #length 2 interactions

#quick fix to get all the length 3 interactions
for inter1 in ppi_nxG.edges_iter("VDR"):
	if inter1[1] == "STAT3":
		print inter1 #0 intermediates
	if inter1[1]==inter1[0]:
		continue
	for inter2 in ppi_nxG.edges_iter(inter1[1]):
		if inter2[1]=="STAT3":
			print inter1[0],inter2[0],inter2[1] #1 intermediates
		for inter3 in ppi_nxG.edges_iter(inter2[1]):
			if inter3[1]=="STAT3":
				print inter1[0],inter2[0],inter3[0],inter3[1]

