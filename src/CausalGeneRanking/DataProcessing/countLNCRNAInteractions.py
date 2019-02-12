"""
	Given all known lncRNAs and a file with all interactions between ncRNA and protein/DNA, how many of these are between lncRNAs and genes?

"""



import sys
import re
import numpy as np

gffInFile = sys.argv[1]
interactionFile = sys.argv[2]


lncRNAs = dict()
with open(gffInFile) as inF:
	
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

lncRNAsWithInteractions = dict()
interactions = 0
with open(interactionFile, 'rb') as inF:
	lineCount = 0
	for line in inF:
		if lineCount < 1:
			lineCount += 1
			continue
		
		splitLine = line.split("\t")
		if splitLine[4] in lncRNAs:
			lncRNAsWithInteractions[splitLine[4]] = 0
			interactions += 1
		
print lncRNAsWithInteractions
print len(lncRNAsWithInteractions)
print interactions