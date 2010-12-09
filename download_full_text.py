import urllib2,os,time,sys,psycopg2
from appscript import *

def articles_from_journal(journal):
	q="""SELECT pmid,medrecord,title,da from "Publication" where journal =%s ORDER BY da"""
	conn=psycopg2.connect("dbname=th17 user=hayssam password=mjnfqto")
	cur=conn.cursor()
	#article with same PMID already here?
	cur.execute(q)
	for r in cur.fetchall():
		pmid,medrecord,title,da=r
		medrecord=eval(medrecord)
		if "AID" in medrecord:
			print da,pmid,title,medrecord["PMID"],medrecord["AID"]
		else:
			print da,pmid,title,"NO AID"

def pmid_of_articles_without_ft_from_journal(journal):
	q="""SELECT pmid FROM "Publication" where journal =%s AND (fulltext IS NULL OR length(fulltext)<=100) ORDER BY da DESC LIMIT 100;"""
	conn=psycopg2.connect("dbname=th17 user=hayssam password=mjnfqto")
	cur=conn.cursor()
	#article with same PMID already here?
	cur.execute(q,(journal,))
	rec=cur.fetchall()
	return [x[0] for x in rec]
		
	
	
def get_doi(pmid):
	rec=get_medline_full_record(pmid)
	doi=None
	if 'AID' not in rec:
		print "No AID for",pmid
		return
	aids=rec["AID"]
	for aid in aids:
		if "doi" in aid:
			doi=aid[:-6]
			break
	return doi




def get_medline_full_record(pmid):
	q="""SELECT medrecord FROM "Publication" where pmid=%s;"""
	conn=psycopg2.connect("dbname=th17 user=hayssam password=mjnfqto")
	cur=conn.cursor()
	#article with same PMID already here?
	cur.execute(q,(pmid,))
	rec=cur.fetchall()[0][0]
	return eval(rec)
	
def store_full_text(pmid,ft):
	q="""UPDATE "Publication" SET fulltext=%s where pmid=%s"""
	conn=psycopg2.connect("dbname=th17 user=hayssam password=mjnfqto")
	cur=conn.cursor()
	cur.execute(q,(ft,pmid))
	conn.commit()
	print "Stored for",pmid,len(ft),"chars"
	sys.stdout.flush()

def print_full_text(pmid):
	q="""SELECT fulltext from "Publication" where pmid=%s"""
	conn=psycopg2.connect("dbname=th17 user=hayssam password=mjnfqto")
	cur=conn.cursor()
	cur.execute(q,(pmid,))
	body=cur.fetchall()[0][0]
	print body
	
def download_process_and_store(pmid,url,delete_before=[],real_start=""):
	f=os.popen("/Volumes/Gauss/bin/lynx -dump \"%s\""%url)
	body=f.read().decode("utf-8",'replace')
	if len(delete_before)>0:
		for marker in delete_before:
			if marker in body:
				start=body.index(marker)+len(marker)
				body=body[start:]
				break
	if real_start!="":
		start=body.index(real_start)+len(real_start)
		related_articles,body=body[:start],body[start:]
		body+=related_articles

	store_full_text(pmid,body)
	
def j_immunol_get_html(pmid):
	url="http://www.jimmunol.org/cgi/pmidlookup?view=long&pmid=%d"%pmid
	f=os.popen("/Volumes/Gauss/bin/lynx -dump \"%s\""%url)
	body=f.read().decode("utf-8",'replace')
	store_full_text(pmid,body)
	sys.stdout.flush()

def get_html_trough_pmc(pmid):
	r=get_medline_full_record(pmid)
	if "PMC" not in r:
		print "No pmc for",pmid
		return
	url="http://www.ncbi.nlm.nih.gov/pmc/articles/%s/?tool=pubmed"%r["PMC"]
	download_process_and_store(pmid,url)
	
sa=app("Safari")
def immunity_get_html(pmid,withSafari=False):
	#get the related pii
	rec=get_medline_full_record(pmid)
	if 'AID' not in rec:
		print "No AID for med",pmid
		return
	aids=rec['AID']
	pii=None
	for aid in aids:
		if 'pii' in aid:
			pii=aid[:-6]
			break
	if pii==None:
		print "No pii for",pmid
		return
	if withSafari:
		url="""http://linkinghub.elsevier.com.gate2.inist.fr/retrieve/pii/%s"""%pii
		t=sa.windows[1].make(new=k.tab)
		sa.windows[1].current_tab.set(t)
		t.URL.set(url)
		time.sleep(2)
		while True:
			try:
				src=t.source()
				if src!=None and src[-8:]=="</html>\n":
					break
			except:
				time.sleep(0.5)
		print "finished loading"
		src=src.encode("utf-8")
		print "Loaded through gate2",len(src),"chars"
		f=open("temp.html","w")
		f.write(src)
		f.close()
		url="temp.html"
	else:
		if pmid in [11970879]:
			url="http://www.sciencedirect.com/science?_ob=ArticleURL&_udi=B6WSP-45N7M44-7&_user=718781&_coverDate=04%2F30%2F2002&_rdoc=1&_fmt=high&_orig=search&_origin=search&_sort=d&_docanchor=&view=c&_acct=C000040138&_version=1&_urlVersion=0&_userid=718781&md5=604e40691e4124fef6fa82ed8c86f7a6&searchtype=a"
		else:
			url="""http://linkinghub.elsevier.com/retrieve/pii/%s"""%pii
	f=os.popen("/Volumes/Gauss/bin/lynx -dump \"%s\""%url)
	body=f.read().decode("utf-8",'replace')
	#move the related articles part from the beg of the body to the end of the body
	marker="Permissions & Reprints\n"
	real_start=body.index(marker)+len(marker)
	related_articles,body=body[:real_start],body[real_start:]
	body+=related_articles
	store_full_text(pmid,body)
	
def eur_j_immunol_get_html(pmid):
	doi=get_doi(pmid)
	if doi==None:
		print "No doi for",pmid
		return
	doi=doi.replace("&#60;","%3C")
	doi=doi.replace("&#62;","%3E")
	url="http://onlinelibrary.wiley.com/doi/%s/full"%doi
	markers=["How to Cite\n","* [48]Permission Request (engl.)","* [33]Podcasts\n"]
	download_process_and_store(pmid,url,delete_before=markers)

def j_exp_med_get_html(pmid):
	url="""http://www.jem.org/cgi/pmidlookup?view=long&pmid=%s"""%pmid
	markers=["* Article\n","* Brief Definitive Report\n"]
	download_process_and_store(pmid,url,delete_before=markers)
	
def blood_get_html(pmid):
	url="""http://www.bloodjournal.org/cgi/pmidlookup?view=long&pmid=%s"""%pmid
	markers=["next article arrow"]
	download_process_and_store(pmid,url,delete_before=markers)

def nature_immunology_get_html(pmid):
	doi=get_doi(pmid)
	if doi==None:
		print "No doi for pmid",pmid
		return
	url="""http://dx.doi.org/%s"""%doi
	download_process_and_store(pmid,url)

def pnas_get_html(pmid):
	url="""http://www.pnas.org/cgi/pmidlookup?view=long&pmid=%s"""%pmid
	markers=["* [15]Sign In as Member / Individual\n"]
	download_process_and_store(pmid,url,delete_before=markers)

def ann_rev_imm_get_html(pmid):
	markers=["right-arrow-small-circle.png"]
	url="http://www.annualreviews.org/doi/full/%s"%get_doi(pmid)
	download_process_and_store(pmid,url,delete_before=markers)

def plos_get_html(pmid):
	url="http://dx.plos.org/"+get_doi(pmid)
	download_process_and_store(pmid,url,delete_before=["\nMetrics"])

def get_oxford(pmid):
	url="http://intimm.oxfordjournals.org/cgi/pmidlookup?view=long&pmid=%s"%pmid
	download_process_and_store(pmid,url)

def get_jci(pmid):
	doi=get_doi(pmid)
	if doi==None:
		print "No doi for pmid",pmid
		return
	url="""http://dx.doi.org/%s"""%doi
	download_process_and_store(pmid,url,delete_before=["Need help?"])

def get_science(pmid):
	url="http://www.sciencemag.org/cgi/pmidlookup?view=full&pmid=%s"%pmid
	download_process_and_store(pmid,url)
	
def get_jbc(pmid):
	url="http://www.jbc.org/cgi/pmidlookup?view=long&pmid=%s"%pmid
	download_process_and_store(pmid,url)

def download_full_text():
	pmids_to_get=pmid_of_articles_without_ft_from_journal("Immunology")
	for pmid in pmids_to_get:
		eur_j_immunol_get_html(pmid)
		time.sleep(0.5)
	print "Finished Immunology"
	
	pmids_to_get=pmid_of_articles_without_ft_from_journal("Science")
	for pmid in pmids_to_get:
		get_science(pmid)
		time.sleep(0.5)
	print "Finished Science"

	pmids_to_get=pmid_of_articles_without_ft_from_journal("Curr Opin Immunol")
	for pmid in pmids_to_get:
		immunity_get_html(pmid)
		time.sleep(0.5)
	print "Finished Curr Opin Immunol"
	
	pmids_to_get=pmid_of_articles_without_ft_from_journal("J Immunol")
	for pmid in pmids_to_get:
		j_immunol_get_html(pmid)
		time.sleep(0.5)
	print "Finished J Immunol"
	
	pmids_to_get=pmid_of_articles_without_ft_from_journal("Blood")
	for pmid in pmids_to_get:
		blood_get_html(pmid)
		time.sleep(1)
	print "Finished Blood"

	pmids_to_get=pmid_of_articles_without_ft_from_journal("Nat Immunol")
	for pmid  in pmids_to_get:
		nature_immunology_get_html(pmid)
		time.sleep(1)
	print "Finished Nat Imm"
	# 
	pmids_to_get=pmid_of_articles_without_ft_from_journal("Proc Natl Acad Sci U S A")
	for pmid  in pmids_to_get:
		pnas_get_html(pmid)
		time.sleep(1)
	print "Finished PNAS"


	pmids_to_get=pmid_of_articles_without_ft_from_journal("Annu Rev Immunol")
	for pmid in pmids_to_get:
		ann_rev_imm_get_html(pmid)
		time.sleep(0.5)
	print "Finished Ann Rev Imm"


	pmids_to_get=pmid_of_articles_without_ft_from_journal("Immunol Rev")
	for pmid in pmids_to_get:
		#wiley article, same as eur j immunol
		eur_j_immunol_get_html(pmid)
		time.sleep(0.5)
	print "Finished Immunol Rev"

	pmids_to_get=pmid_of_articles_without_ft_from_journal("Clin Exp Immunol")
	for pmid in pmids_to_get:
		#wiley article, same as eur j immunol
		eur_j_immunol_get_html(pmid)
		time.sleep(0.5)
	print "Finished Clin Exp Immunol"


	pmids_to_get=pmid_of_articles_without_ft_from_journal("Nature")
	for pmid in pmids_to_get:
		#wiley article, same as eur j immunol
		nature_immunology_get_html(pmid)
		time.sleep(0.5)
	print "Finished Nature"


	pmids_to_get=pmid_of_articles_without_ft_from_journal("Arthritis Rheum")
	for pmid in pmids_to_get:
		#wiley article, same as eur j immunol
		eur_j_immunol_get_html(pmid)
		time.sleep(0.5)
	print "Finished Arthritis Rheum"

	pmids_to_get=pmid_of_articles_without_ft_from_journal("PLoS One")
	for pmid in pmids_to_get:
		plos_get_html(pmid)
		time.sleep(0.5)
	print "Finished PLoS One"
	
	pmids_to_get=pmid_of_articles_without_ft_from_journal("Nat Rev Immunol")
	for pmid in pmids_to_get:
		#wiley article, same as eur j immunol
		nature_immunology_get_html(pmid)
		time.sleep(0.5)
	print "Finished Nat Rev Immunol"


	pmids_to_get=pmid_of_articles_without_ft_from_journal("Nat Med")
	for pmid in pmids_to_get:
		#wiley article, same as eur j immunol
		nature_immunology_get_html(pmid)
		time.sleep(0.5)
	print "Finished Nat Med"

	pmids_to_get=pmid_of_articles_without_ft_from_journal("Cell")
	for pmid in pmids_to_get:
		immunity_get_html(pmid)
		time.sleep(0.5)
	print "Finished Cell"

	pmids_to_get=pmid_of_articles_without_ft_from_journal("Int Immunol")
	for pmid in pmids_to_get:
		get_oxford(pmid)
		time.sleep(0.5)
	print "Finished Int Immunol"
	
	pmids_to_get=pmid_of_articles_without_ft_from_journal("J Clin Invest")
	for pmid in pmids_to_get:
		get_jci(pmid)
		time.sleep(0.5)
	print "Finished J Clin Invest"

	pmids_to_get=pmid_of_articles_without_ft_from_journal("J Allergy Clin Immunol")
	for pmid in pmids_to_get:
		immunity_get_html(pmid,withSafari=True)
		time.sleep(0.5)
	print "Finished J Allergy Clin Immunol"
	
	pmids_to_get=pmid_of_articles_without_ft_from_journal("Mol Immunol")
	for pmid in pmids_to_get:
		immunity_get_html(pmid,withSafari=True)
		time.sleep(0.5)
	print "Finished Mol Immunol"

	pmids_to_get=pmid_of_articles_without_ft_from_journal("J Biol Chem")
	for pmid in pmids_to_get:
		get_jbc(pmid)
		time.sleep(0.5)
	print "Finished J Biol Chem"

	
def mark_reviews():
	q="""SELECT pmid,title,medrecord from "Publication" where medrecord LIKE '%Review%' ORDER BY da"""
	qupd="""UPDATE "Publication" SET isreview=%s WHERE pmid=%s"""
	conn=psycopg2.connect("dbname=th17 user=hayssam password=mjnfqto")
	cur=conn.cursor()
	#article with same PMID already here?
	cur.execute(q)
	for r in cur.fetchall():
		pm,title,medrecord=r
		medrecord=eval(medrecord)
		if "PT" in medrecord and "Review" in medrecord["PT"]:
			cur2=conn.cursor()
			cur2.execute(qupd,(True,pm))
		elif "STAT" in medrecord and "Review" in medrecord["STAT"]
	conn.commit()
		