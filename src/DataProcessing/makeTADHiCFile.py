"""
	Make a file with Hi-C interactions for every TAD. Each side of the interaction is considered a separate regulatory element that can be gained or lost.
	
	1. Filter for all interactions that do not take place entirely within a TAD
	2. Write these interactions per TAD to a file.

"""

import sys
import numpy as np

hiCFile = sys.argv[1]
tadFile = sys.argv[2]
outFile = sys.argv[3]

#1. Read the TADs

tads = []
with open(tadFile, 'r') as inF:
	
	for line in inF:
		line = line.strip()
		splitLine = line.split("\t")
		
		tads.append([splitLine[0], int(splitLine[1]), int(splitLine[2])])

	tads = np.array(tads, dtype="object")



#2. Read the Hi-C interactions and link to TADs. 
tadInteractions = dict()
with open(hiCFile, 'r') as inF:
	lineCount = 0
	for line in inF:
		
		line = line.strip()
		splitLine = line.split("\t")
		
		region1 = splitLine[0]
		region2 = splitLine[1]
		
		regionChr = region1.split("_")[0]
		region1Pos = region1.split("_")[1]
		region2Pos = region2.split("_")[1]
		
		#1. Map both regions to a TAD.
		
		tadChrSubset = tads[tads[:,0] == regionChr]
		#Match region 1
		matchingStart = tadChrSubset[:,1] <= int(region1Pos)
		matchingEnd = tadChrSubset[:,2] >= int(region1Pos)
		
		matchingPos = matchingStart * matchingEnd
		if len(tadChrSubset[matchingPos]) < 1:
			continue
		matchingTad1 = tadChrSubset[matchingPos][0] 
		
		#Match region 1
		matchingStart = tadChrSubset[:,1] <= int(region2Pos)
		matchingEnd = tadChrSubset[:,2] >= int(region2Pos)
		
		matchingPos = matchingStart * matchingEnd
		if len(tadChrSubset[matchingPos]) < 1:
			continue
		matchingTad2 = tadChrSubset[matchingPos][0]
		
		#Only consider interactions that are within the same TAD. We assume that interactions always take place within TADs, so otherwise this doesn't fit with our data.
		
		if matchingTad1[0] == matchingTad2[0] and matchingTad1[1] == matchingTad2[1] and matchingTad1[2] == matchingTad2[2]:
			tadStr = matchingTad1[0] + "_" + str(matchingTad1[1]) + "_" + str(matchingTad1[2])
			if tadStr not in tadInteractions:
				tadInteractions[tadStr] = []
			tadInteractions[tadStr].append(region1Pos)
			tadInteractions[tadStr].append(region2Pos)
		lineCount += 1

#Write the Hi-c interaction indices to a file grouped by TADs
with open(outFile, 'w') as outF:
	for tad in tadInteractions:
		
		line = tad
		line += "\t"
		for interaction in tadInteractions[tad]:
			line += str(interaction)
			line += ","
		line += "\n"
		outF.write(line)
	


