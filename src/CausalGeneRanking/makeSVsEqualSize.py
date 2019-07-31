"""
	For every somatic SV, get the germline SV of the same type with the most similar size.
	Then shrink or enlarge that SV to the same size (on both ends)
	Output these SVs to a new file. 

"""

import sys
import numpy as np
from inputParser import InputParser
import math
from random import randrange


somaticSVs = np.load(sys.argv[1])
germlineSVs = np.load(sys.argv[2])

sizedGermlineSVs = []

#Make subsets
deletions = []
for germlineSV in germlineSVs:
	typeMatch = False
	if germlineSV[8].svType == "deletion":
		deletions.append(germlineSV)
deletions = np.array(deletions,dtype="object")

duplications = []
for germlineSV in germlineSVs:
	typeMatch = False
	if germlineSV[8].svType == "duplication":
		duplications.append(germlineSV)
duplications = np.array(duplications,dtype="object")

inversions = []
for germlineSV in germlineSVs:
	typeMatch = False
	if germlineSV[8].svType == "inversion":
		inversions.append(germlineSV)
inversions = np.array(inversions,dtype="object")

for somaticSV in somaticSVs:
	
	somaticSize = somaticSV[5] - somaticSV[1]
	#print "somatic SV", somaticSize
	
	currentSmallestSizeDiff = float("inf")
	currentGermlineSize = float("inf")
	currentBestSV = None
	
	
	#First subsample the right type
	matchingGermlines = []
	
	if somaticSV[8].svType == "del":
		matchingGermlines = deletions
	if somaticSV[8].svType == "tandem_dup":
		matchingGermlines = duplications
	if somaticSV[8].svType == "invers":
		matchingGermlines = inversions
	
	if len(matchingGermlines) == 0:
		continue
	
	randomSV = randrange(matchingGermlines.shape[0])
	#print "random ind: ", randomSV
	ind = 0
	for germlineSV in matchingGermlines:
		
		#Determine if this SV has a size more similar than any other SV of this type
		size = germlineSV[5] - germlineSV[1]
		
		#Select the size of a random SV, not the nearest!
		# if ind == randomSV:
		# 	currentBestSV = germlineSV
		# 	currentGermlineSize = np.abs(somaticSize - size)
		# 	#print "germline SV: ", currentGermlineSize
		# 	break
		# 
		# ind += 1
		
		if np.abs(somaticSize - size) < currentSmallestSizeDiff:
			currentSmallestSizeDiff = np.abs(somaticSize - size)
			currentBestSV = germlineSV
			currentGermlineSize = size
		
	#From the SV with the most similar size, make the size completely equal
	
	if currentGermlineSize == float("inf"): #SV type has no match with any germline, skip this
		continue
	
	if somaticSize > currentGermlineSize: #if the somatic SV is larger, enlarge the germline SV on both sides
		
		sizeDifference = somaticSize - currentGermlineSize
		halfSizeDifference = np.round(sizeDifference / 2.0)
		
		newSV = [somaticSV[0], int(somaticSV[1]+halfSizeDifference), int(somaticSV[2]+halfSizeDifference), somaticSV[8].o1,
				 somaticSV[3], int(somaticSV[4]-halfSizeDifference), int(somaticSV[5]-halfSizeDifference), somaticSV[8].o2,
				 'brca_tcga_equalSize', somaticSV[7], somaticSV[8].svType, somaticSV[6]]
		
		if newSV[4] < 0:
			newSV[4] = 0
		if newSV[5] < 0:
			newSV[5] = 0
		
		sizedGermlineSVs.append(newSV)
	
	if somaticSize <= currentGermlineSize: #if the somatic SV is larger, enlarge the germline SV on both sides
		sizeDifference = currentGermlineSize - somaticSize
		halfSizeDifference = np.round(sizeDifference / 2.0)
		
		newSV = [somaticSV[0], int(somaticSV[1]-halfSizeDifference), int(somaticSV[2]-halfSizeDifference), somaticSV[8].o1,
				 somaticSV[3], int(somaticSV[4]+halfSizeDifference), int(somaticSV[5]+halfSizeDifference), somaticSV[8].o2,
				 'brca_tcga_equalSize', somaticSV[7], somaticSV[8].svType, somaticSV[6]]
		
		if newSV[1] < 0:
			newSV[1] = 0
		if newSV[2] < 0:
			newSV[2] = 0
		
		sizedGermlineSVs.append(newSV)
		
header = 'chr1	s1	e1	o1	chr2	s2	e2	o2	source	sample_name	sv_type	cancer_type'		
		
np.savetxt('Output/somaticSVs_equalSized.txt', np.array(sizedGermlineSVs, dtype="object"), fmt="%s", delimiter='\t', header=header)
	
	

	
	
	
	

