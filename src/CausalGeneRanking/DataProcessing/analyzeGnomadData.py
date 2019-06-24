"""
	Analysis of gnomad data. Type distribution, length, etc. 

"""
import sys
import numpy as np

gnomadData = np.loadtxt(sys.argv[1], dtype="object")

print gnomadData
header = dict()
with open(sys.argv[1], 'r') as inF:
	lineCount = 0
	for line in inF:
		if lineCount < 1:
			line = line.strip()
			splitLine = line.split("\t")
			for ind in range(0, len(splitLine)):
				val = splitLine[ind]
				header[val] = ind
			break
print header	

#How many of each type are in the data?

types = dict()

for svType in gnomadData[:,header['SVTYPE']]:
	if svType not in types:
		types[svType] = 0
	types[svType] += 1

import matplotlib.pyplot as plt

# plt.bar(types.keys(), types.values())
# plt.show()

#What is the length distribution of the SVs?

# svLens = gnomadData[:,header['SVLEN']]
# 
# plt.hist(svLens)
# plt.show()

#How many SVs do not PASS?

filterTypes = dict()

for filterType in gnomadData[:,header['FILTER']]:
	if filterType not in filterTypes:
		filterTypes[filterType] = 0
	filterTypes[filterType] += 1

print filterTypes
# plt.bar(filterTypes.keys(), filterTypes.values())
# plt.show()

#Filter and see how many SVs remain
filteredSVs = []
for sv in gnomadData:
	if sv[header['FILTER']] == "PASS":
		#if sv[header['SVTYPE']] == 'DEL' or sv[header['SVTYPE']] == 'DUP' or sv[header['SVTYPE']] == 'INV' or sv[header['SVTYPE']] == 'CTX':
		if sv[header['SVTYPE']] == 'DEL' or sv[header['SVTYPE']] == 'DUP' or sv[header['SVTYPE']] == 'INV':
			filteredSVs.append(sv)

filteredSVs = np.array(filteredSVs, dtype="object")
print filteredSVs.shape
#exit()

#How many TADs are overlapped by these germline SVs?

tadData = np.loadtxt(sys.argv[2], dtype='object')
overlappedTADCounts = dict()
for sv in filteredSVs:
	
	tadChrSubset = tadData[tadData[:,0] == 'chr' + sv[0]]
	
	startMatches = float(sv[1]) <= tadChrSubset[:,2].astype('float')
	endMatches = float(sv[2]) >= tadChrSubset[:,1].astype('float')
	
	tadMatches = tadChrSubset[startMatches * endMatches]
	# 
	# if tadMatches.shape[0] > 100:
	# 	print sv
	# 	exit()

	if tadMatches.shape[0] not in overlappedTADCounts:
		overlappedTADCounts[tadMatches.shape[0]] = 0
	overlappedTADCounts[tadMatches.shape[0]] += 1	

print overlappedTADCounts	
#plt.hist(overlappedTADCounts)
#plt.show()

	
	


#Is there any overlap with the somatic SVs? 

