"""
	For every somatic SV, get the germline SV of the same type with the most similar size.
	Then shrink or enlarge that SV to the same size (on both ends)
	Output these SVs to a new file. 

"""

import sys
import numpy as np
from inputParser import InputParser
import math


somaticSVs = np.load(sys.argv[1])
germlineSVs = np.load(sys.argv[2])

sizedGermlineSVs = []

for somaticSV in somaticSVs:
	
	somaticSize = somaticSV[5] - somaticSV[1]
	
	currentSmallestSizeDiff = float("inf")
	currentGermlineSize = float("inf")
	currentBestSV = None
	for germlineSV in germlineSVs:
		#Match the right type
		typeMatch = False
		
		if somaticSV[8].svType == "del" and germlineSV[8].svType == "deletion":
			typeMatch = True
		if somaticSV[8].svType == "tandem_dup" and germlineSV[8].svType == "duplication":
			typeMatch = True
		if somaticSV[8].svType == "invers" and germlineSV[8].svType == "inversion":
			typeMatch = True	

		if typeMatch == False:
			continue
		
		#Determine if this SV has a size more similar than any other SV of this type
		size = germlineSV[5] - germlineSV[1]
		
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
		
		newSV = [currentBestSV[0], int(currentBestSV[1]-halfSizeDifference), int(currentBestSV[2]-halfSizeDifference), currentBestSV[8].o1,
				 currentBestSV[3], int(currentBestSV[4]+halfSizeDifference), int(currentBestSV[5]+halfSizeDifference), currentBestSV[8].o2,
				 'germline_dgv_equalSize', currentBestSV[7], currentBestSV[8].svType, currentBestSV[6]]
		
		sizedGermlineSVs.append(newSV)
	
	if somaticSize <= currentGermlineSize: #if the somatic SV is larger, enlarge the germline SV on both sides
		sizeDifference = currentGermlineSize - somaticSize
		halfSizeDifference = np.round(sizeDifference / 2.0)
		
		newSV = [currentBestSV[0], int(currentBestSV[1]+halfSizeDifference), int(currentBestSV[2]+halfSizeDifference), currentBestSV[8].o1,
				 currentBestSV[3], int(currentBestSV[4]-halfSizeDifference), int(currentBestSV[5]-halfSizeDifference), currentBestSV[8].o2,
				 'germline_dgv_equalSize', currentBestSV[7], currentBestSV[8].svType, currentBestSV[6]]
		
		sizedGermlineSVs.append(newSV)
		
header = 'chr1	s1	e1	o1	chr2	s2	e2	o2	source	sample_name	sv_type	cancer_type'		
		
np.savetxt('Output/germlineSVs_equalSized.txt', np.array(sizedGermlineSVs, dtype="object"), fmt="%s", delimiter='\t', header=header)
	
	

	
	
	
	

