"""
	For every TAD boundary, plot where the nearest SVs are in the somatic case and the germline case.
	To show this, limit the counts to the next TAD boundary. 

"""

import sys
import numpy as np
import matplotlib.pyplot as plt
from inputParser import InputParser


somaticSVs = np.load(sys.argv[1])

filteredSomaticSVs = [] #first filter the SVs to remove the translocations, there are not relevant here
for somaticSV in somaticSVs:
	if somaticSV[0] == somaticSV[3]:
		filteredSomaticSVs.append(somaticSV)
filteredSomaticSVs = np.array(filteredSomaticSVs, dtype="object")

tads = InputParser().getTADsFromFile(sys.argv[2])

#for every TAD boundary, count until the next TAD boundary on each side how many SV breakpoints are there.
coordinates = dict()
for tad in tads:
	
	
	#start with tad start. Find the nearest TAD boundary on the left side.
	#The nearest tad boundary on the right will be the end of this TAD.
	nearestLeftFromStart = []
	nearestLeftDist = float("inf")
	nearestRightFromEnd = []
	nearestRightDist = float("inf")
	for tad2 in tads:
		if tad2[3] == tad[3]:
			continue
		
		if tad[0] != tad2[0]: #don't compare TADs on different chromosomes
			continue
		
		startDistance = tad[1] - tad2[2]
		if startDistance < 0:
			continue #This TAD is not on the left
		if startDistance < nearestLeftDist:
			nearestLeftDist = startDistance
			nearestLeftFromStart = tad2
		
		endDistance = tad2[1] - tad[2]
		if endDistance < 0:
			continue #This TAD is not on the right
		if endDistance < nearestRightDist:
			nearestRightDist = endDistance
			nearestRightFromEnd = tad2
	
	#if either value is empty, this is the farmost left or right TAD, so look at all SVs smaller or larger than that.
	#thus the left dist can be set  to 0, and the right can stay infinite,assuming that SVs don't go outside chromosomes
	if len(nearestLeftFromStart) == 0:
		nearestLeftDist = 0
	
	maxWindow = 1000
	
	
	### TAD START ###
	
	#Get all SVs from the nearest left until the TAD start
	
	svChrMatches = filteredSomaticSVs[filteredSomaticSVs[:,0] == tad[0]]
	
	startMatches = svChrMatches[(svChrMatches[:,1] >= nearestLeftDist) * (svChrMatches[:,1] <= tad[1])]
	endMatches = svChrMatches[(svChrMatches[:,5] >= nearestLeftDist) * (svChrMatches[:,5] <= tad[1])]
	
	#Get all the positions at which these breakpoints are and use these as x coordinates
	for match in startMatches:
		bp = match[1]
		#The distance is relative to the TAD start
		distance = bp - tad[1]
		
		if distance < -maxWindow: #exclude positions too far away.
			continue
		if distance not in coordinates:
			coordinates[distance] = 0
		coordinates[distance] += 1
	
	for match in endMatches:
		bp = match[5]
		#The distance is relative to the TAD start
		distance = bp - tad[1]
		if distance < -maxWindow: #exclude positions too far away.	
			continue
		if distance not in coordinates:
			coordinates[distance] = 0
		coordinates[distance] += 1
	#Repeat but then for the right side of the right TAD boundary
	
	svChrMatches = filteredSomaticSVs[filteredSomaticSVs[:,0] == tad[0]]
	
	startMatches = svChrMatches[(svChrMatches[:,1] <= tad[2]) * (svChrMatches[:,1] >= tad[1])]
	endMatches = svChrMatches[(svChrMatches[:,5] <= tad[2]) * (svChrMatches[:,5] >= tad[1])]
	
	#Get all the positions at which these breakpoints are and use these as x coordinates
	for match in startMatches:
		bp = match[1]
		#The distance is relative to the TAD start
		distance = bp - tad[1]
		if distance > maxWindow: #exclude positions too far away. 
			continue
		if distance not in coordinates:
			coordinates[distance] = 0
		coordinates[distance] += 1 
	
	for match in endMatches:
		bp = match[5]
		#The distance is relative to the TAD start
		distance = bp - tad[1]
		if distance > maxWindow: #exclude positions too far away. 
			continue
		if distance not in coordinates:
			coordinates[distance] = 0
		coordinates[distance] += 1

	### TAD END ###
	
	#Get all SVs from the tad start until the tad end
	#This is the same as above, so we can repeat that step
	
	svChrMatches = filteredSomaticSVs[filteredSomaticSVs[:,0] == tad[0]]
	
	startMatches = svChrMatches[(svChrMatches[:,1] >= tad[1]) * (svChrMatches[:,1] <= tad[2])]
	endMatches = svChrMatches[(svChrMatches[:,5] >= tad[1]) * (svChrMatches[:,5] <= tad[2])]
	
	#Get all the positions at which these breakpoints are and use these as x coordinates
	for match in startMatches:
		bp = match[1]
		#The distance is relative to the TAD end
		distance = bp - tad[2]
		if distance < -maxWindow: #exclude positions too far away. 
			continue
		if distance not in coordinates:
			coordinates[distance] = 0
		coordinates[distance] += 1
	
	for match in endMatches:
		bp = match[5]
		#The distance is relative to the TAD start
		distance = bp - tad[2]
		if distance < -maxWindow: #exclude positions too far away. 
			continue
		if distance not in coordinates:
			coordinates[distance] = 0
		coordinates[distance] += 1
	
	#Repeat but then for the right side of the right TAD boundary
	
	svChrMatches = filteredSomaticSVs[filteredSomaticSVs[:,0] == tad[0]]
	
	startMatches = svChrMatches[(svChrMatches[:,1] <= nearestRightDist) * (svChrMatches[:,1] >= tad[2])]
	endMatches = svChrMatches[(svChrMatches[:,5] <= nearestRightDist) * (svChrMatches[:,5] >= tad[2])]
	
	#Get all the positions at which these breakpoints are and use these as x coordinates
	for match in startMatches:
		bp = match[1]
		#The distance is relative to the TAD start
		distance = bp - tad[2]
		if distance > maxWindow: #exclude positions too far away.
			
			continue
		if distance not in coordinates:
			coordinates[distance] = 0
		coordinates[distance] += 1
	
	for match in endMatches:
		bp = match[5]
		#The distance is relative to the TAD start
		distance = bp - tad[2]
		if distance > maxWindow: #exclude positions too far away. 
			continue
		if distance not in coordinates:
			coordinates[distance] = 0
		coordinates[distance] += 1

#Bin the values for visualization
print max(coordinates.keys())
print min(coordinates.keys())

print coordinates

binSize = 5
binnedCoordinates = dict()
sortedCoordinates = np.sort(coordinates.keys())
currentBinStart = sortedCoordinates[0]
binnedCoordinates = dict()
currentBin = int(round(-((max(coordinates.keys()) - min(coordinates.keys())) / float(binSize)) / 2)) #oof
for coordinate in sortedCoordinates:
	
	if coordinate < (currentBinStart + binSize):
		if currentBin not in binnedCoordinates:
			binnedCoordinates[currentBin] = 0
		binnedCoordinates[currentBin] += 1
	else:
		currentBinStart += binSize
		currentBin += 1

print binnedCoordinates
print "plotting"		
plt.bar(binnedCoordinates.keys(), binnedCoordinates.values())
plt.show()