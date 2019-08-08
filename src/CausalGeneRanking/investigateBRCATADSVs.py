"""
	Some relevant information for the disruption of the TAD/TAD boundary next to BRCA1.

"""

#Read all SVs

from __future__ import absolute_import
from __future__ import print_function
from inputParser import InputParser
from tad import TAD
from gene import Gene
import settings

import numpy as np
import matplotlib.pyplot as plt

svFile = settings.files['svFile']
svData = InputParser().getSVsFromFile(svFile, "all")

#Define the relevant TAD
brcaTADRight = TAD("chr17", 41440000, 41480000)
brcaGene = Gene("BRCA1", "chr17", 41197696, 41276113)

#Get all the SVs that overlap with this TAD

overlappingSVs = []
overlappingSVTypes = dict()

tadOverlappingSVs = []
tadOverlappingSVTypes = dict()

allPatients = dict()
allPatientsSVs = dict()

patientsWithSV = dict()

for sv in svData:
	if sv[7] not in allPatients:
		allPatients[sv[7]] = 0
		allPatientsSVs[sv[7]] = []
	
	#start of sv before TAD end, end of sv after TAD start
	if sv[0] != brcaTADRight.chromosome:
		continue

	if sv[1] < brcaTADRight.end and sv[5] > brcaTADRight.start:
		overlappingSVs.append(sv)
		if sv[8].svType not in overlappingSVTypes:
			overlappingSVTypes[sv[8].svType] = 0
		overlappingSVTypes[sv[8].svType] += 1
		
	#Filtering out the SVs that also overlap with BRCA1:
	if sv[1] < brcaTADRight.end and sv[5] > brcaTADRight.start and sv[1] > brcaGene.start:
		tadOverlappingSVs.append(sv)
		if sv[8].svType not in tadOverlappingSVTypes:
			tadOverlappingSVTypes[sv[8].svType] = 0
		tadOverlappingSVTypes[sv[8].svType] += 1
		allPatients[sv[7]] += 1
		allPatientsSVs[sv[7]].append([sv[0], sv[1], sv[5], sv[8].svType])
		patientsWithSV[sv[7]] = 0
		
# 		
# print overlappingSVTypes
# print tadOverlappingSVTypes
# 
# print allPatients
# print np.unique(np.array(tadOverlappingSVs)[:,7])

print(np.unique(list(patientsWithSV.keys())))
print(len(list(allPatients.keys())))



#What are the types of the overlapping SVs? 

