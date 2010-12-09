#!/usr/bin/env python
# encoding: utf-8
import copy

# project="Th17 differentiation"
# termFile=project+"/articles_th17.txt"
# pubmed_files=[f for f in os.listdir("./"+project) if f.startswith("entrez ")]

stop_words=open("stop_words.txt").readlines()
stop_words=[x.strip() for x in stop_words]

occ_number={}
exemple_sentences={}
n_grams_occurences={}

def build_n_grams(n,table):
	n_grams=[]
	for i in range(0,len(table)-n+1):
		n_grams.append(tuple(table[i:i+n]))
	return n_grams
	
# toRemove=re.compile("|".join(["of","the","in","for","with","by","this","well","as"]))

def tokenize_sentence(s):
	s=s.replace("(","")
	s=s.replace(")","")
	s=s.replace(",","")

	splitted=[]
	recognized=allTermsRe.findall(s)
	for parts in allTermsRe.split(s):
		if (parts in recognized):
			parts=parts.lower()
			if "_"+parts in allTerms:
				parts="_"+parts
			splitted.append(allTerms[parts].name)
		else:
			for w in parts.split():
				if  w.lower() in stop_words:
					continue
				w=w.lower()
				w=w.strip(",.")
				splitted.append(w.lower())
	return splitted
			
			
def compute_occurences():
	global occ_number,n_grams_occurences
	occ_number={}
	n_grams_occurences={}
	
	for s in allSentences:
		#we get rid of ()
		s=s.string
		splitted=[]
		tokens=tokenize_sentence(s)
		for tok in tokens:
			if tok not in occ_number:
				occ_number[tok]=0
			occ_number[tok]+=1
			if tok not in exemple_sentences:
				exemple_sentences[tok]=[]
			exemple_sentences[tok].append(s)
		for i in range(2,6):
			grams=build_n_grams(i,tokens)
			for g in grams:
				g=tuple(g)
				if g not in n_grams_occurences:
					n_grams_occurences[g]=0
				if g not in exemple_sentences:
					exemple_sentences[g]=[]
				n_grams_occurences[g]+=1
				exemple_sentences[g].append(s)
# print occ_number


#we first build the frequencies of terms related to all publications
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

project="Th17 differentiation"
termFile="Th17 differentiation human/th17_human_terms.txt"
pubmed_files=[f for f in os.listdir("./"+project) if f.startswith("entrez ")]
parse_file(termFile)

compute_occurences()
occ_number_general=copy.copy(occ_number)
n_grams_occurences_general=copy.copy(n_grams_occurences)
exemple_sentences_general=copy.copy(exemple_sentences)
n_sents_general=len(allSentences)
print "General term frequencies computed"


allPublications={}
allSentences=[]
allPredicates=[]
pubmedIdTopub={}
nxG=None

project="Th17 differentiation human"
termFile=project+"/th17_human_terms.txt"
pubmed_files=[f for f in os.listdir("./"+project) if f.startswith("entrez ")]
pubmed_files.append("pubmed_result.txt")
parse_file(termFile)
occ_number={}
exemple_sentences={}
n_grams_occurences={}
compute_occurences()
n_sents=len(allSentences)
print "Human specific frequencies computed"

diff_freq={}
for k,v in occ_number.items():
	if k not in diff_freq:
		diff_freq[k]=0
	if k not in occ_number_general:
		diff_freq[k]=1.0*v/n_sents
	else:
		diff_freq[k]=1.0*occ_number[k]/n_sents-1.0*occ_number_general[k]/n_sents_general

#sort the dicts by freqs
occ_number=occ_number.items()
occ_number.sort(key=lambda x:x[1],reverse=True)

occ_number_general=occ_number_general.items()
occ_number_general.sort(key=lambda x:x[1],reverse=True)

n_grams_occurences=[x for x in n_grams_occurences.items() if x[1]>3]
n_grams_occurences.sort(key=lambda x:x[1],reverse=True)

n_grams_occurences_general=[x for x in n_grams_occurences_general.items() if x[1]>3]
n_grams_occurences_general.sort(key=lambda x:x[1],reverse=True)

diff_freq=diff_freq.items()
diff_freq.sort(key=lambda x:x[1],reverse=True)

