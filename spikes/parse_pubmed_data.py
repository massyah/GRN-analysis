#!/usr/bin/env python
# encoding: utf-8
import re

#get all abstracts
abstract_re=re.compile("AB(.*?)AD",re.DOTALL)
pmid_id_re=re.compile("PMID\- ([0-9]*)")
first_author=re.compile("AD.*?FAU(.*?)AU",re.DOTALL)
title_re=re.compile("TI(.*?)AB",re.DOTALL)

txt=file("TGF-β - IL23r.txt","r").read()
# txt=txtTemplate

# print txt[:2000]
# for m in abstract_re.findall(txt):
# 	print m[4:20]

for m in title_re.findall(txt):
	print m