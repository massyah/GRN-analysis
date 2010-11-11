#!/usr/bin/env python
# encoding: utf-8

"""
Generic combinatorial functions
* Generate all subsets
"""
def gensubs(n):
	yield []
	for i in range(1, 2 ** n):
		yield filter(lambda j: i & 2 ** j, range(n))
		
def allsets(n): 
	"""Read more: http://adorio-research.org/wordpress/?p=6791#ixzz13S0NoKiS"""
	yield []
	S = [[i] for i in range(n)] #initial. 
	for s in S:
		yield s 
	N = S 
	for k in range(n): 
		if N: 
			S = N 
			N = [] 
			for e in S: 
				for i in range(e[-1]+1,n): 
					f = e[:] 
					f.append(i) 
					yield f
					N.append(f) 
def mysubsets(l):
	for i in range(0,2**len(l)):
		items=[]
		idx=0
		while i >0:
			y=i%2
			i/=2
			if y:
				items.append(l[idx])
			idx+=1
		#pad
		# binrep+="0"*(n-len(binrep))
		yield items


# print len([i for i in gensubs(10)]) #0.09 secs, 1024 subsets
# print len([i for i in gensubs(18)]) #2.56 secs, 262144 subsets
# gensubs(32) #fails range(2**n) is not possible
# subs= [i for i in gensubs(3)] #0.09 secs, 8 subsets
# print len(subs),subs


subs=[x for x in allsets(20)] #0.94 secs, 262144 subsets
print len(subs)
# subs=[x for x in allsets(5)] #0.94 secs, 262144 subsets
# print len(subs)

# l=range(20)
# subs=[x for x in mysubsets(l)]
# print len(subs)
# for s in subs:
# 	print ",".join(s)