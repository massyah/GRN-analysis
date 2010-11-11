#we describe a small part of activations in the th17 model and output that to ginSim

project="Th17 differentiation"
termFile=project+"/articles_th17.txt"
pubmed_files=[]


global allPublications, allSentences, allTerms, allTermsRe, allPredicates, uid, sentencesUid, evidencesUid, pubmedIdTopub, nxG,cmd
allPublications={}
allSentences=[]
allTerms={}
allTermsRe=None
allPredicates=[]
uid=0
sentencesUid=0
evidencesUid=0
pubmedIdTopub={}
nxG=None
parse_file()

cmd=generate_variable_command("""ifng,ifngr,stat1,socs1,il-4,il-4r,stat6,gata3,t-bet,il-12,il-12r,stat4,il-18,il-18r,irak,ifnb,ifnbr,il21,stat3,tgfb,ahr,
il-1b, il-12, il-4, il-6, il-5, il-18, il-2, il-21, il-23,il-17,il-23r,il-27
il-2r,il-6r,il-21r,il-27r
smad3, nfat,foxp3,stat5,jak3,foxp3prom,runx1,socs3,jak,stat1,jak1,jak2
il17a,il17f,rorgt,il21,il21r,il-1b,il-1br,il22
cd4p,tcr,nfkb,ap1,ets1,ifngprom,rora
myd88,
traf6,
jnk,
ikk,
jun,
fos,
jak2-jak1,
tyk2-jak2,
jak1-jak3,
ibp,irf4,smad4

""")
exec(compile(cmd,"","exec"))


template="""<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE gxl SYSTEM "http://gin.univ-mrs.fr/GINsim/GINML_2_1.dtd">
<gxl xmlns:xlink="http://www.w3.org/1999/xlink">
	<graph id="default_name" class="regulatory" nodeorder="%(nodeList)s">
	%(nodeDefs)s
	%(edgesDefs)s
	</graph>
	</gxl>"""
nodeTemplate="""<node id="%(id)s" basevalue="%(basal)i" maxvalue="%(max)i">
"""
nodeParams="""<parameter idActiveInteractions="%(edgesId)s" val="%(focal)i"/>"""
edgeTemplate="""<edge id="%(id)s" from="%(src)s" to="%(tgt)s" minvalue="%(activationLevel)i" sign="%(sign)s">\n</edge>\n"""

#signaling axis from the ILR
SignalingAxis([il_2,il_2r,jak1_jak3,stat5])
BindsRepress(stat5,il17a,15)
Activates(il17a,il_17,36)
BindsRepress(foxp3,il17a,36)

build_nx_graph(onlyPositive=False)
# print nxG.edges(data=True)

nodes=nxG.nodes(data=True)
nodeList=" ".join([x[0] for x in nodes])
nodeDefs=""
for n in nodes:
	nodeDefs+=nodeTemplate%{"id":n[0],"basal":0,"max":1}
	influencing=nxG.predecessors(n[0])
	if len(influencing):
		for i in influencing:
			nodeType=type(predicates_with_uid(nxG[i][n[0]]["uid"])[0])
			if nodeType in [Downregulates,Inhibits,Inactivates,Sequestrates,BindsRepress]:
				focalVal=0
			else:
				focalVal=1
			nodeDefs+=nodeParams%{
			"edgesId":"_".join([i,n[0]]),
			"focal":focalVal
			}
	nodeDefs+="\n</node>\n"
	
edgeDefs=""
def add_activation(e):
	return edgeTemplate%{
	"id":"_".join(e[:2]),
	"src":e[0],
	"tgt":e[1],
	"activationLevel":1,
	"sign":"positive"
	}
def add_inhibition(e):
	print "adding inhib"
	return edgeTemplate%{
	"id":"_".join(e[:2]),
	"src":e[0],
	"tgt":e[1],
	"activationLevel":1,
	"sign":"negative"
	}
	

for e in nxG.edges(data=True):
	print e
	p=predicates_with_uid(e[2]["uid"])[0]
	if type(p) in [Activates,Upregulates,Binds]:
		f=add_activation
	elif type(p) in [Downregulates,Inhibits,Inactivates,Sequestrates,BindsRepress]:
		f=add_inhibition
	else:
		f=add_activation
	edgeDefs+=f(e)

print edgeDefs
f=open("generatedModel.ginml","w")
f.write(template%{
	"nodeList":nodeList,
	"nodeDefs":nodeDefs,
	"edgesDefs":edgeDefs
	})
f.close()