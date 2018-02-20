#!/usr/bin/env python

######### Documentation #########

__author__ = "Marleen Nieboer"
__credits__ = []
__maintainer__ = "Marleen Nieboer"
__email__ = "m.m.nieboer@umcutrecht.nl"
__status__ = "Development"

"""
	This is a very basic script for investigating the quality of the True Positve SV dataset. It can check for:
		- Recurrence within the true positive set
		- Overlap with the true negative set

"""

### Imports ###

import numpy as np
import sys
import matplotlib.pyplot as plt

### Code ###


#1. Check how many other SVs overlap with each SV in the true positive set (recurrence)

recurrence = dict()
overlapWindow = 0 #windw of additional bp overlap to check for
exactOverlap = False #If true, we look for SVs within the defined window of bp of other SVs. Otherwise, any overlap is allowed. 

# Read the data into a format with only the relevant position information
svData = []
inFile = sys.argv[1]
	
with open(inFile, "r") as f:
	
	lineCount = 0
	header = []
	for line in f:
		line = line.strip()
		splitLine = line.split("\t")
		if lineCount < 1:
			header = splitLine
			lineCount += 1
			continue
		
		#First, we check if the header contains SV-related values or not. If not, it is an SNV/SNP.
		
		#Obtain required data for each data type. Search which column has this information
		chr1Index = header.index("chr1")
		s1Index = header.index("s1")
		e1Index = header.index("e1")
		chr2Index = header.index("chr2")
		s2Index = header.index("s2")
		e2Index = header.index("e2")
		
		chr1 = splitLine[chr1Index]
		chr2 = splitLine[chr2Index]
		s1 = int(splitLine[s1Index])
		e1 = int(splitLine[e1Index])
		
		if splitLine[s2Index] != 'NaN':
			s2 = int(splitLine[s2Index])
		else:
			s2 = s1
		if splitLine[e2Index] != 'NaN':
			e2 = int(splitLine[e2Index])
		else:
			e2 = e1
		
		
		svData.append([chr1,s1,e1,chr2,s2,e2])
		
#First sort the SVs to make sure that subsetting is as efficient as possible
svDataNp = np.array(svData, dtype='object')
	
chr1col = svDataNp[:,0]
chr2col = svDataNp[:,3]
	
sortedInd = np.lexsort((chr2col, chr1col)) #sort first by column 1, then by column 2. This works, but it is lexographical, so chromosome 11 comes before chromosome 2. For this purpose it is ok, since we only look for this
												#chromosome in other files to subset these, so the order in which we do that does not really matter. 	
sortedSvData = svDataNp[sortedInd]

#Load the true negatives as well if we wish to compare to this dataset
if len(sys.argv) > 2:
	trueNegatives = sys.argv[2]
	trueNegativesData = []
		
	with open(trueNegatives, "r") as f:
		
		lineCount = 0
		header = []
		for line in f:
			line = line.strip()
			splitLine = line.split("\t")
			if lineCount < 1:
				header = splitLine
				lineCount += 1
				continue
			
			#First, we check if the header contains SV-related values or not. If not, it is an SNV/SNP.
			
			#Obtain required data for each data type. Search which column has this information
			chr1Index = header.index("chr1")
			s1Index = header.index("s1")
			e1Index = header.index("e1")
			chr2Index = header.index("chr2")
			s2Index = header.index("s2")
			e2Index = header.index("e2")
			
			chr1 = splitLine[chr1Index]
			chr2 = splitLine[chr2Index]
			s1 = int(splitLine[s1Index])
			e1 = int(splitLine[e1Index])
			
			if splitLine[s2Index] != 'NaN':
				s2 = int(splitLine[s2Index])
			else:
				s2 = s1
			if splitLine[e2Index] != 'NaN':
				e2 = int(splitLine[e2Index])
			else:
				e2 = e1
			
			trueNegativesData.append([chr1,s1,e1,chr2,s2,e2])
	comparisonDataset = np.array(trueNegativesData, dtype='object')
	
else:
	comparisonDataset = sortedSvData #here we can easily switch out for another dataset, the script will do the same


#Now perform overlap operations
#For each SV, check if there are other SVs that are within the defined overlap window

previousChr1 = None
previousChr2 = None
overlapDistribution = dict()
for svNum in range(0, len(sortedSvData)):
	
	#First make a subset of chromosomes to search through, reduce the search space
	lineList = sortedSvData[svNum,:]

	#We should check the chromosome of the previous line.
	#This would work nicely if the input file has been properly sorted! We probably need a pre-sorting to make this work. 
	if str(lineList[0]) != previousChr1 and str(lineList[3]) != previousChr2:
		
		#Find the two subsets that match on both chromosomes.
		matchingChr1Ind = comparisonDataset[:,0] == lineList[0]
		matchingChr2Ind = comparisonDataset[:,0] == lineList[3] #The gene list entry is here the same since the genes are only on 1 chromosome. 
		
		#The SV is only the same if both chromosomes are the same, so we do not need to look at individual chromosome subsets. 
		matchingChrInd = matchingChr1Ind * matchingChr2Ind
		chrSubset = comparisonDataset[np.where(matchingChrInd)]
		
		#Make sure to update the previous chromosome when it changes
		previousChr1 = str(lineList[0])
		previousChr2 = str(lineList[3])
	
	if np.size(chrSubset) < 1:
		continue #no need to compute the distance, there are no genes on these chromosomes
	
	if lineList[0] == lineList[3]: #If the chromosomes are the same
		start = lineList[1]#only do the conversion once. Why is this not remembered from the previous loop?
		end = lineList[5]
	else: #if the chromosomes are not the same, then we only use the positions of chromosome 1. 
		start = lineList[1]
		end = lineList[2]
		
	if exactOverlap == True:
		startOverlap1 = start - overlapWindow <= chrSubset[:,1]
		startOverlap2 = start + overlapWindow >= chrSubset[:,1]
		endOverlap1 = end - overlapWindow <= chrSubset[:,5]
		endOverlap2 = end + overlapWindow >= chrSubset[:,5]
		
		matches = startOverlap1 * startOverlap2 * endOverlap1 * endOverlap2
	else:
		startOverlap = start - overlapWindow <= chrSubset[:,5]
		endOverlap = end + overlapWindow >= chrSubset[:,1]
	
		#Now find where both of these are true (multiply because both conditions need to be true)
		matches = startOverlap * endOverlap
	
	matchIndChr = np.where(matches == True)[0] #This creates an array in array, we are only interested in the first that contains the matching indices
	
	matchingSVs = chrSubset[matchIndChr, :]
	
	#This is the number of SVs that each SV overlaps with within the same dataset.
	if len(matchingSVs) > 0:
		overlapDistribution[svNum] = len(matchingSVs)
	
	
print overlapDistribution

confidentlyRecurrent = 0
for sv in overlapDistribution.keys():
	
	if overlapDistribution[sv] > 1:
		confidentlyRecurrent += 1

print "Number of SVs that are confidently recurring more than once within a window of : ", overlapWindow, " bp: ", confidentlyRecurrent
#Visualize the overlap	

plt.bar(overlapDistribution.keys(), overlapDistribution.values(), 50)
plt.xlabel("SVs")
plt.ylabel("Frequency")
axes = plt.gca()
plt.show()
	
	
	
	




		
			
			
		