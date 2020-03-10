#check how many SVs have gains/losses, and use this to report how many SVs disrupt TADs.

import numpy as np
import sys

rules = np.loadtxt(sys.argv[1], dtype='object')

uniqueSVs = dict()
svTypes = dict()
for sv in rules[:,0]:

	splitSV = sv.split('_')
	svStr = '_'.join(splitSV[1:])
	
	svType = splitSV[12]
	
	if svStr not in uniqueSVs:
		uniqueSVs[svStr] = 0
		if svType not in svTypes:
			svTypes[svType] = 0
		svTypes[svType] += 1
	
	
print(len(uniqueSVs))
print(svTypes)


gl = np.loadtxt(sys.argv[2], dtype='object')

uniqueSVs = 0
svTypes = dict()
lineCount = 0
for sv in gl:
	
	if lineCount < 1:
		lineCount += 1
		continue

	uniqueSVs += 1

	svType = sv[10]
	if svType not in svTypes:
		svTypes[svType] = 0
	svTypes[svType] += 1
	
	
print(uniqueSVs)
print(svTypes)