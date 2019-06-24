# Get the p-values per gene, assign score of 1 to signficant and 0 to non-significant

import numpy as np
import sys

pValues = np.loadtxt(sys.argv[1], dtype="object")

binPValues = []
for pValue in pValues[:,1]:
	if float(pValue) == 0 or float(pValue) < 0.05:
		binPValues.append(1)
	else:
		binPValues.append(0)

newPValues = np.empty(pValues.shape, dtype="object")

newPValues[:,0] = pValues[:,0]
newPValues[:,1] = binPValues

np.savetxt(sys.argv[2], newPValues, delimiter='\t', fmt='%s')

