import re

txt=open("1generatedModel.ginml","r").read()
nodePos_re=re.compile("<node id=\"(.+?)\".*?<nodevisualsetting>(.*?)</nodevisualsetting>",re.DOTALL)
edgePos_re=re.compile("<edge .*? from=\"(.+?)\" to=\"(.+?)\".*?<edgevisualsetting>(.*?)</edgevisualsetting>",re.DOTALL)

pos={}
edgePos={}
for m in nodePos_re.finditer(txt):pos[m.groups()[0]]=m.groups()[1].strip()

for m in edgePos_re.finditer(txt): edgePos[m.groups()[0],m.groups()[1]]=m.groups()[2].strip()

print pos
print edgePos