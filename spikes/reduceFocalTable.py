g="Wg6"
focal=focal_table(g)

#we take the rows != 1
active=[x for x in focal.items() if x[1]!=0]

influence=influencing_genes(g)
#we associate to each var in the table, the value it can take to activate the tgt
newTable={}
for gene in influence:
	possibleValues=set(range(domain(gene)+1))
	for a in active:
		d=dict(a[0])
		if d[gene] in possibleValues:
			possibleValues.remove(d[gene])
	for v in possibleValues:
		row=dict(zip(influence,"*"*len(influence)))
		row[gene]=v
		newTable[tuple(row.items())]=0

for k,v in active:
	newTable[k]=v
basal=tuple(zip(influence,[0]*len(influence)))
newTable[basal]=0
return newTable