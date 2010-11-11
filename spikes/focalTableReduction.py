ft={(('IL4', 0), ('SOCS1', 0)): 0,
 (('IL4', 0), ('SOCS1', 1)): 0,
 (('IL4', 1), ('SOCS1', 0)): 1,
 (('IL4', 1), ('SOCS1', 1)): 0}

# def domain(g):
# 	if g=="STAT1":
# 		return 2
# 	else:
# 		return 1


def reduce_focal_table(ft):
	genes=set()
	for row in ft.items():
		for col in row[0]:
			genes.add(col[0])
		
	inputs=set([x[0] for x in ft.items()])
	newRows={}
	for g in genes:
		for i in range(domain(g)+1):
			val=(g,i)
			focals=set([x[1] for x in ft.items() if val in x[0]])
			if len(focals)==1:
				v=focals.pop()
				print "Reduce on",val,"to",v
				inpToRemove=[x[0] for x in ft.items() if val in x[0]]
				inputs=inputs-set(inpToRemove)
				newRows[(val)]=v
	for i in inputs:
		newRows[i]=ft[i]
	return newRows
			
print reduce_focal_table(ft)
