


import sys

sys.path.insert(0, '../')

import numpy as np
import re
import matplotlib.pyplot as plt
from inputParser import InputParser

# svFile = sys.argv[1]
# 
# variantsList = []
# 
# with open(svFile, 'rb') as f:
# 	
# 	lineCount = 0
# 	header = []
# 	for line in f:
# 		line = line.strip()
# 		splitLine = line.split("\t")
# 		
# 		#First extract the header and store it in the dictionary to remove dependency on the order of columns in the file
# 		if lineCount < 1:
# 
# 			header = splitLine
# 			lineCount += 1
# 			continue
# 		
# 		#Now extract the chromosome, start and end (there are multiple, 2 for an SV)
# 		chr1Index = header.index("chr1")
# 		s1Index = header.index("s1")
# 		e1Index = header.index("e1")
# 		o1Index = header.index("o1")
# 
# 		chr2Index = header.index("chr2")
# 		s2Index = header.index("s2")
# 		e2Index = header.index("e2")
# 		o2Index = header.index("o2")
# 
# 		cancerTypeIndex = header.index("cancer_type")
# 		sampleNameIndex = header.index("sample_name")
# 		
# 		cancerType = splitLine[cancerTypeIndex]
# 		sampleName = splitLine[sampleNameIndex]
# 		
# 		# #Dirty workaround to make sure that the cancer type names are the same, we only focus on 1 type for this intial run
# 		if cancerType == "breast/gastric":
# 			cancerType = "breast"
# 			
# 		#Skip anything that is not breast cancer for now. From here is the easiest way, saves time in processing as well
# 		# if cancerType != settings.general['cancerType']:
# 		# 	continue
# 		
# 		svTypeIndex = header.index("sv_type")
# 		svType = splitLine[svTypeIndex]
# 		
# 		# if typeFilter != "all":
# 		# 	#Check if the SV type matches deletions
# 		# 	match = re.search("deletion", svType, re.IGNORECASE)
# 		# 	if match is None: #only focus on deletions for now
# 		# 		continue
# 		# 
# 		
# 		#if cancerType not in uniqueCancerTypes:
# 		#	uniqueCancerTypes.append(cancerType)
# 
# 		#If the coordinates are missing on the second chromosome, we use the coordinates of the first chromosome unless chr 1 and chr 2 are different.
# 		if splitLine[chr1Index] == splitLine[chr2Index]:
# 			if splitLine[s2Index] == 'NaN':
# 				splitLine[s2Index] = int(splitLine[s1Index])
# 				
# 			if splitLine[e2Index] == 'NaN':
# 				splitLine[e2Index] = int(splitLine[e1Index])
# 		else:
# 			if splitLine[chr2Index] == 'NaN':
# 				continue # This line does not have correct chromosome 2 information (should we be skipping it?)
# 
# 		s1 = int(splitLine[s1Index])
# 		e1 = int(splitLine[e1Index])
# 		s2 = int(splitLine[s2Index])
# 		e2 = int(splitLine[e2Index])
# 		chr2 = splitLine[chr2Index]
# 		
# 		chr1 = splitLine[chr1Index]
# 		o1 = splitLine[o1Index]
# 		o2 = splitLine[o2Index]
# 		
# 		#Make sure to switch the positions here as well
# 		#Some positions are swapped
# 		if int(e2) < int(e1):
# 			tmpE1 = e1
# 			e1 = e2
# 			e2 = tmpE1
# 			tmpS1 = s1
# 			s1 = s2
# 			s2 = tmpS1
# 		
# 		#Sometimes only the end is swapped.
# 		if int(e2) < int(s2):
# 			tmpS2 = s2
# 			s2 = e2
# 			e2 = tmpS2
# 			
# 		if int(e1) < int(s1):
# 			tmpS1 = s1
# 			s1 = e1
# 			e1 = tmpS1
# 		
# 		
# 		#svObject = SV('chr' + chr1, s1, e1, o1, 'chr' + chr2, s2, e2, o2, sampleName, cancerType, svType)
# 		#chr 1, start, end, chr2, start2, end2
# 		variantsList.append(['chr' + chr1, s1, e1, 'chr' + chr2, s2, e2, cancerType, sampleName, svType])
# 
# regions = np.array(variantsList, dtype='object')
regions = InputParser().getSVsFromFile(sys.argv[1], "all")


#Plot the SV size per type
deletions = []
inversions = []
duplications = []
for sv in regions:
	if sv[8].svType == "deletion":
		deletions.append(sv)

	if sv[8].svType == "invers":
		inversions.append(sv)
		
	if sv[8].svType == "tandem_dup":
		duplications.append(sv)

deletions = np.array(deletions, dtype="object")
inversions = np.array(inversions, dtype="object")
duplications = np.array(duplications, dtype="object")

sizes = []
for deletion in inversions:
	size = deletion[5] - deletion[1]
	sizes.append(size)

plt.hist(sizes)
plt.show()


# 
# #Determine how many TADs are overlapped by each SV
# tadFile = sys.argv[2]
# tads = []
# with open(tadFile, 'r') as tadF:
# 	
# 	for line in tadF:
# 		line = line.strip()
# 		splitLine = line.split("\t")
# 		
# 		tads.append([splitLine[0], int(splitLine[1]), int(splitLine[2])])
# 	
# tads = np.array(tads, dtype="object")
# 
# #Go through deletions, inv and dups, to skip translocations
# overlappedTadCounts = dict()
# sizes = []
# for sv in deletions:
# 	
# 	tadSubset = tads[tads[:,0] == sv[0]]
# 	
# 	startMatches = sv[1] < tadSubset[:,1]
# 	endMatches = sv[5] > tadSubset[:,2]
# 	
# 	allMatches = startMatches * endMatches
# 	overlappedTads = tadSubset[allMatches]
# 	
# 	
# 	
# 	if overlappedTads.shape[0] > 2 or sv[5] - sv[1] > 30000000: #skip large SV on chry
# 		continue
# 	
# 	if sv[5] - sv[1] > 20000000:
# 		print sv
# 		print overlappedTads
# 		print overlappedTads.shape
# 		
# 
# 	sizes.append(sv[5] - sv[1])
# 
# 	if len(overlappedTads) not in overlappedTadCounts:
# 		overlappedTadCounts[len(overlappedTads)] = 0
# 	overlappedTadCounts[len(overlappedTads)] += 1
# 
# print len(sizes)
# print np.min(sizes)
# print np.max(sizes)
# print np.mean(sizes)
# exit()
# 	
# 
# print overlappedTadCounts
# plt.bar(overlappedTadCounts.keys(), overlappedTadCounts.values())
# plt.show()
# exit()


#Plot in bins of predefined sizes (< 1kb, )

#Deletions
sizes = []
sizes = dict()
sizes[10] = 0
sizes[10**2] = 0
sizes[10**3] = 0
sizes[10**4] = 0
sizes[10**5] = 0
sizes[10**6] = 0
sizes[10**7] = 0
sizes[10**8] = 0
sizes[10**9] = 0

# 
# for deletion in deletions:
# 	size = deletion[5] - deletion[1]
# 	#sizes.append(size)
# 	if size > 10 and size < 10**2:
# 		sizes[10] += 1
# 	if size > 10**2 and size < 10**3:
# 		sizes[10**2] += 1
# 	if size > 10**3 and size < 10**4:
# 		sizes[10**3] += 1
# 	if size > 10**4 and size < 10**5:
# 		sizes[10**4] += 1
# 	if size > 10**5 and size < 10**6:
# 		sizes[10**5] += 1
# 	if size > 10**6 and size < 10**7:
# 		sizes[10**6] += 1
# 	if size > 10**7 and size < 10**8:
# 		sizes[10**7] += 1
# 	if size > 10**8 and size < 10**9:
# 		sizes[10**8] += 1
# 	if size > 10**9 and size < 10**10:
# 		sizes[10**9] += 1


for deletion in deletions:
	size = deletion[5] - deletion[1]
	#sizes.append(size)
	if size > 10 and size < 10**2:
		sizes[10] += 1
	if size > 10**2 and size < 10**3:
		sizes[10**2] += 1
	if size > 10**3 and size < 10**4:
		sizes[10**3] += 1
	if size > 10**4 and size < 10**5:
		sizes[10**4] += 1
	if size > 10**5 and size < 10**6:
		sizes[10**5] += 1
	if size > 10**6 and size < 10**7:
		sizes[10**6] += 1
	if size > 10**7 and size < 10**8:
		sizes[10**7] += 1
	if size > 10**8 and size < 10**9:
		sizes[10**8] += 1
	if size > 10**9 and size < 10**10:
		sizes[10**9] += 1
# 
print sizes
exit()
# 
# print np.max(sizes)
# print np.min(sizes)
# print np.mean(sizes)


#Plot histogram
plt.bar(sizes.keys(), sizes.values())
plt.show()
# 
# (n, bins, patches) = plt.hist(sizes, bins=[10, 10**2, 10**3, 10**4, 10**5, 10**6, 10**7, 10**8, 10**9])
# print n
# plt.show()
#Inversions


#Duplications

