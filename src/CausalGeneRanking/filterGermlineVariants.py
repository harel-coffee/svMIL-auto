"""

Subsample germline variants to the same size as HMF variants. 

"""
import sys
import numpy as np


variantCount = 73293

germlineFile = sys.argv[1]

germlineData = []
with open(germlineFile, 'r') as inF:
	lineCount = 0
	for line in inF:
		
		if lineCount < 1:
			
			lineCount += 1
			continue

		splitLine = line.split('\t')
		
		if splitLine[5] not in ['deletion', 'tandem duplication', 'inversion']:
			continue
		
		if len(splitLine[1]) > 2: #skip weird chromosomes
			continue
		
		chr1 = 'chr' + splitLine[1]
		
		if splitLine[5] == 'deletion':
			svType = 'DEL'
		elif splitLine[5] == 'tandem duplication':
			svType = 'DUP'
		else:
			svType = 'INV'
		
		#use old format for parsing with tool
		#chr1, s1, e1, o1, chr2, s2, e2, o2, source, sample anme, sv type, cancer type
		germlineData.append([chr1, int(splitLine[2]), int(splitLine[2]), '.', chr1,
							 int(splitLine[3]), int(splitLine[3]), '.', 'DGV', 'na',
							 svType, 'germline'])

germlineData = np.array(germlineData, dtype='object')

#random subsample

randomIndices = np.random.randint(germlineData.shape[0], size=variantCount)

randomIndices = np.sort(randomIndices)

randomGermlineData = germlineData[randomIndices,:]

header = 'chr1\ts1\te1\to1\tchr2\ts2\te2\to2\tsource\tsample_name\tsv_type\tcancer_type'

with open(sys.argv[2], 'w') as outF:
	
	outF.write(header)
	outF.write('\n')
	
	for row in range(0, randomGermlineData.shape[0]):
		rowLine = str(randomGermlineData[row][0])
		
		for col in range(1, randomGermlineData.shape[1]):
			
			rowLine += '\t'
			rowLine += str(randomGermlineData[row][col])
			
		outF.write(rowLine)
		outF.write('\n')
	