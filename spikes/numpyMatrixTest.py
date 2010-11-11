import numpy,random
def test_col_assignment():
	a=numpy.zeros((1000,1000))
	a[3,5]=11
	a[5,6]=22
	a[:,5]=55
	a[:,6]=1
	print a

#test speed vs python arrays

def test_1d_array_npy():
	N=65
	a=numpy.zeros(65)
	for i in range(len(a)):
		a[i]=1
	return a
def test_1d_array_npy_row():
	N=65
	a=numpy.zeros(65)
	a[0:]=1
	return a
def test_1d_array_py():
	N=65
	a=[0]*65
	for i in range(len(a)):
		a[i]=1
	return a

def test_1d_array(): #slower for npy
	for i in range(20000):
		test_1d_array_py()

def mod_array_py(a):
	for i in range(len(a)):
		a[i]=1
def mod_array_npy_row(a):
	a[0:]=1

def test_mod_array(): #slower for npy
	a=[0]*65
	an=numpy.array(a)
	for i in range(200000):
		mod_array_py(a)
	
def test_sum_vectors_npy(a,b):
	return a+b
def test_sum_vectors_py(a,b):
	c=[0]*len(a)
	for i in range(len(a)):
		c[i]=a[i]+b[i]
	return c

def test_sum():#slower for npy, except on big arrays (size > 10), then way faster
	a=[1,2,3]*100
	b=[1,2,3]*100
	an=numpy.array(a)
	bn=numpy.array(b)
	for i in range(200000):
		test_sum_vectors_npy(an,bn)
		# test_sum_vectors_py(a,b)


Any=-3 #-1*the max possible value of any state of the model -1 
geneRange=range(4)
def merge_states(d1,d2):
	global ncall
	d1=list(d1)
	for i in geneRange:
		if d1[i]!=d2[i]:
			if d1[i]==Any:
				d1[i]=d2[i]
			elif d2[i]==Any:
				continue
			else:
				return None
	return tuple(d1)
# 
def merge_states_npy(d1,d2):
	res=d1+d2+3
	if numpy.sometrue(res>2):
		return None
	else:
		return res

s1=[Any,Any,0,2]*50
s2=[Any,2,Any,Any]*50
s3=[Any,Any,1,2]*50
np=map(numpy.array,[s1,s2,s3])

def test_merge():
	for i in range(20000):
		merge_states(s1,s2)
		merge_states(s1,s3)
def test_merge_npy():
	for i in range(20000):
		merge_states_npy(np[0],np[1])
		merge_states_npy(np[0],np[2])

#reprent the state_merge with a * operator?, see also the other operators, like sum
#see the ix_() function for fuse_table
#ndenumerate()
#ndindex()
#unique() faster that set?


#swith to md5.new(X).digest for storing rows without duplicates?,http://www.mail-archive.com/numpy-discussion@scipy.org/msg16491.html