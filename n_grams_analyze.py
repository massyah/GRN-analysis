#companion to n_grams, assume it has been runned before running this

def exemples(k):
	"""k can be a single term or a n-gram"""
	sents=exemple_sentences[k]
	sents=[Sentence(None,s) for s in sents]
	for s in sents:
		print s.highlight()

#for the words having different frequencies, we print those who are not in the /usr/share/dict

common_english=open("/usr/share/dict/words").readlines()
common_english=[x.strip() for x in common_english]

def domain_specific_words_freq_diff():
	for k,v in diff_freq[:500]:
		if type(k)==str and k not in common_english:
			print k,v

