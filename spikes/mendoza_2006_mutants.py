"""
We automate the verification of the WT and mutant phenotypes by using the stables states
We read the model from the ginml file
"""
from cProfile import run as profilerun
import numpy



from IPython.Debugger import Tracer; debug_here = Tracer()

import ginSimParser_003 as gp

print gp.reducedFocalParamsValues


import fuse_tables_008 as ft

#to get the phenotypes st states profiles: 
#for i in range(len(ft.tables[0])): print ft.tables[0][i][1:]

wt=({},[
(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
(0, 0, 2, 1, 0, 0, 0, 0, 1, 1, 2, 0, 0, 0, 0, 0, 0),
(0, 0, 1, 1, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0),
(0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0)
])
il_12_p=({"IL12":[1,1]},[
(0, 0, 2, 1, 0, 0, 0, 1, 1, 1, 2, 1, 0, 0, 1, 0, 0),
(0, 0, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 1, 0, 0),
(0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0)
])

ifng_0=({"IFNg":[0,0]},[
(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
(0, 0, 2, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0),
(0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0),
(0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0)
]
)


ifngr_0=({"IFNgR":[0,0]},[
(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
(0, 0, 2, 0, 0, 0, 0, 0, 1, 0, 2, 0, 0, 0, 0, 0, 0),
(0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0),
(0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0)
])

il_18_p=({"IL18":[1,1]},[
(0, 1, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0),
(1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1),
(1, 0, 2, 1, 0, 1, 0, 0, 1, 1, 2, 0, 0, 0, 0, 0, 1),
(1, 0, 1, 1, 0, 1, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 1)
])

il_12_18_p=({"IL12":[1,1],"IL18":[1,1]},[
(1, 0, 2, 1, 0, 1, 0, 1, 1, 1, 2, 1, 0, 0, 1, 0, 1),
(1, 0, 1, 1, 0, 1, 0, 1, 1, 1, 2, 1, 0, 0, 1, 0, 1),
(0, 1, 0, 0, 1, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0)
])
gata3_p=({"GATA3":[1,1]},[
(0, 0, 2, 1, 1, 0, 0, 0, 1, 1, 2, 0, 0, 0, 0, 0, 0),
(0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0),
(0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0)
])


mutationsPhenotypes=[wt,il_12_p,ifng_0,ifngr_0,il_18_p,il_12_18_p,gata3_p]

for m in mutationsPhenotypes:
	tables=ft.compute_sstate(gp.reducedFocalParamsValues,m[0])
	print m[0],
	valid=True
	for ssState in tables[0]:
		valid = valid & (ssState[1:] in m[1])
	print valid
	print "-"*8