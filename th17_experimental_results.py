import os,itertools,unicodedata,sys

project="Th17 differentiation"
termFile=project+"/articles_th17.txt"
pubmed_files=[f for f in os.listdir("./"+project) if f.startswith("entrez ")]
pubmed_files=[]

model_name="th17.ginml"

from IPython.Debugger import Tracer; debug_here = Tracer()


import omg_interface as omg

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
all_edges_id=[]
parse_file()

cmd=generate_variable_command("""ifng,ifngr,stat1,socs1,il-4,il-4r,stat6,gata3,t-bet,il-12,il-12r,stat4,il-18,il-18r,irak,ifnb,ifnbr,il21,stat3,tgfb,ahr,
il-1b, il-12, il-4, il-6, il-5, il-18, il-2, il-21, il-23,il-17,il-23r,il-27
il-2r,il-6r,il-21r,il-27r
smad3, nfat,foxp3,stat5,jak3,foxp3prom,runx1,socs3,jak,stat1,jak1,jak2
il17a,il17f,rorgt,il21,il21r,il-1b,il-1br,il22,il23r
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
import os,itertools,unicodedata,sys

project="Th17 differentiation"
termFile=project+"/articles_th17.txt"
pubmed_files=[f for f in os.listdir("./"+project) if f.startswith("entrez ")]
pubmed_files=[]

model_name="th17.ginml"

from IPython.Debugger import Tracer; debug_here = Tracer()


import omg_interface as omg

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
all_edges_id=[]
parse_file()

cmd=generate_variable_command("""ifng,ifngr,stat1,socs1,il-4,il-4r,stat6,gata3,t-bet,il-12,il-12r,stat4,il-18,il-18r,irak,ifnb,ifnbr,il21,stat3,tgfb,ahr,
il-1b, il-12, il-4, il-6, il-5, il-18, il-2, il-21, il-23,il-17,il-23r,il-27
il-2r,il-6r,il-21r,il-27r
smad3, nfat,foxp3,stat5,jak3,foxp3prom,runx1,socs3,jak,stat1,jak1,jak2
il17a,il17f,rorgt,il21,il21r,il-1b,il-1br,il22,il23r
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

mine=Publication("My findings","hayssam",2010)
ra=allTerms["_ra"]
il17_secretion=allTerms["il-17 secretion"]
il21_secretion=allTerms["il-21 secretion"]