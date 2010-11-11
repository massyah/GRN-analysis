#we describe a small part of activations in the th17 model and output that to ginSim
import unicodedata

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

SignalingAxis([il_2,il_2r,jak1_jak3,stat5],oshea2009)
SignalingAxis([il_2,il_2r,jak1_jak3,stat3],oshea2009)
SignalingAxis([il_2,il_2r,jak1_jak3,stat1],oshea2009)
SignalingAxis([il_6,il_6r,jak2_jak1,stat3],7)
SignalingAxis([il_21,il_21r,jak1_jak3,stat3],7)
SignalingAxis([il_21,il_21r,jak1_jak3,stat1],7)
SignalingAxis([il_21,il_21r,jak1_jak3,stat5],7)
SignalingAxis([il_23,il_23r,tyk2_jak2,stat3],7)
SignalingAxis([il_23,il_23r,tyk2_jak2,stat4],14)
SignalingAxis([il_1b,il_1br,myd88,irak,traf6],27)



Activates([smad3,nfat],foxp3,[43,18157133])
Activates(tgfb,smad3) #very old result
Binds(stat3,il17a)
Binds(stat3,il17f,8)
Binds(stat5,foxp3prom)
Activates(tgfb,foxp3) #from PMID 14676299
Upregulates(t_bet,il_23r,1)
Activates(tcr,nfat,2)
Binds(nfat,il17a,3)
Activates(tcr,nfkb,4)
Activates(tcr,ap1,4)
Complexes(ap1,nfat,5)
Binds(stat4,ifng)
Binds(stat3,il21,9)
Upregulates(stat3,il_23r,10)

Downregulates(socs3,stat3,11) #also see <38>, socs3 negreg il-23-mediated stat3 signaling, does it holds for other il-n mediated stat3 activation?

Activates(ifng,stat1,12)
Upregulates([il_27,ifng],t_bet,13)
d=BindsRepress(stat5,il17a,15) #seem validated in pmid 17363300
d.hypothesis=True

Activates(allTerms["il-2 secretion"],il_2)
Inhibits(ets1,allTerms["il-2 secretion"],19)

Binds(stat5,foxp3prom)
Inhibits(foxp3,rorgt)
Inhibits(rorgt,foxp3)
Upregulates(il21r,il_21r)
Upregulates(foxp3prom,foxp3)
Activates(stat3,il21r,16)

Downregulates(stat5,rorgt,[17,20696842])

Binds([ets1,t_bet],ifngprom,18)
Upregulates(ifngprom,ifng,direct=True)

Activates(tgfb,smad3,20)
Activates([tgfb,il_6],rorgt,21)
Activates([tgfb,il_2],foxp3,21)

Upregulates([tgfb,il_6],il_23r,22)
Inhibits(stat3,foxp3,23)
Downregulates(tgfb,il_23r,24)
Downregulates(tgfb,t_bet,25)
Downregulates(tgfb,allTerms["il-2 secretion"],25)
Activates([il_1b,il_6],rorgt)
Activates(traf6,jnk,28)
Activates(traf6,ikk,28)
Phosphorylates(jnk,jun,29)
Activates([fos,jun],ap1,29)
Activates(ikk,nfkb,30)
Binds(rorgt,il17a,31)
Upregulates(stat3,rorgt,32)
Activates([stat3,rorgt],il17_secretion,33)
Activates(il17a,il17_secretion)
Activates(il17f,il17_secretion)
RequiredToActivate([rora,rorgt],il17_secretion,34)
Downregulates(ra,rorgt,35)
Upregulates(ra,foxp3,35)
BindsRepress(foxp3,il17a,36)
Downregulates(foxp3,il17_secretion,37)
Downregulates(foxp3,il21_secretion,37)
Upregulates(il21,il21_secretion)
Sequestrates(ibp,irf4,39)
Upregulates([irf4,stat3],rorgt,40)
Binds([runx1,rorgt],il17a,41)
Sequestrates(foxp3,runx1,17377532) #yields runx1 that can't dimerize with rorgt
Activates(tcr,ibp,15249594)
Upregulates([il_21,rorgt,stat3],il_23r,44)
Inhibits(smad3,il_2,1388)

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
	"id":"_".join(e[:2]).replace(" ","_"),
	"src":e[0],
	"tgt":e[1],
	"activationLevel":1,
	"sign":"positive"
	}
def add_inhibition(e):
	return edgeTemplate%{
	"id":"_".join(e[:2]).replace(" ","_"),
	"src":e[0],
	"tgt":e[1],
	"activationLevel":1,
	"sign":"negative"
	}
	

for e in nxG.edges(data=True):
	# print e
	p=predicates_with_uid(e[2]["uid"])[0]
	if type(p) in [Activates,Upregulates,Binds]:
		f=add_activation
	elif type(p) in [Downregulates,Inhibits,Inactivates,Sequestrates,BindsRepress]:
		f=add_inhibition
	else:
		f=add_activation
	edgeDefs+=f(e)


#HACK bug with ANDs, 
model
# print edgeDefs
f=open("generatedModel.ginml","w")
model=template%{
	"nodeList":nodeList,
	"nodeDefs":nodeDefs,
	"edgesDefs":edgeDefs
	}
model=model.replace("AND ","AND_")
#remove extended unicode chars to play well with gin sim

model=unicodedata.normalize('NFKC', model).encode('ascii','ignore')
f.write(model)
f.close()