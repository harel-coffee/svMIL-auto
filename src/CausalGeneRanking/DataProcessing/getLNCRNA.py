"""
	The goal of this script is to filter the lncRNAs from the refseq file with all genomic information.
"""

import sys
import re
import numpy as np

inFile = sys.argv[1]
outFile = sys.argv[2]


lncRNAs = dict()
with open(inFile) as inF:
	
	for line in inF:
		
		if re.search("#", line) == None:
			line = line.strip()
			splitLine = line.split("\t")
			
			if splitLine[2] == "lnc_RNA":
				
				info = splitLine[8]
				splitInfo = info.split(";")
				fieldDict = dict()
				for field in splitInfo:
					splitField = field.split("=")
					
					fieldDict[splitField[0]] = splitField[1]
				
				
					
				if fieldDict["gene"] not in lncRNAs:
					lncRNAs[fieldDict["gene"]] = []
				#Complicated chromosome encoding
				splitLocusId = splitLine[0].split(".")
				
				splitSplitLocusId = splitLocusId[0].split("_")
				
				#Now replace all 0s in the second part
				chromId = splitSplitLocusId[1].replace("0", "")
				
				
				
				
				lncRNAs[fieldDict["gene"]].append([fieldDict["gene"], 'chr' + chromId, splitLine[3], splitLine[4]])
				
			
#Filter for the longest transcript	
filteredLncRNAs = []			
for lncRNA in lncRNAs:
	
	if len(lncRNAs[lncRNA]) > 1: #if there are more genes with the same name, get the longest
		currentLongest = lncRNAs[lncRNA][0]
		currentLength = int(lncRNAs[lncRNA][0][3]) - int(lncRNAs[lncRNA][0][2])
		for transcript in lncRNAs[lncRNA]:
			length = int(transcript[3]) - int(transcript[2])
			if length > currentLength:
				currentLength = length
				currentLongest = transcript
		
		filteredLncRNAs.append(currentLongest)
	else:
		filteredLncRNAs.append(lncRNAs[lncRNA][0])
		
print len(filteredLncRNAs)
print np.array(filteredLncRNAs)

#Write the genes to an output file

with open(outFile, 'w') as outF:
	
	for lncRNA in filteredLncRNAs:
		outF.write(lncRNA[1] + "\t" + lncRNA[2] + "\t" + lncRNA[3] + "\t" + lncRNA[0] + "\n")
		
	
		