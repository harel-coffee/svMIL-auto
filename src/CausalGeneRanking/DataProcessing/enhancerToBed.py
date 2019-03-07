import sys
import re
enhancerFile = sys.argv[1]
outFile = sys.argv[2]

enhancers = []
with open(outFile, 'w') as outF:
	with open(enhancerFile, 'rb') as f:
		
		lineCount = 0
		for line in f:
			if lineCount < 1:
				lineCount += 1
				continue
			
			line = line.strip()
			splitLine = line.split("\t")
			
			interaction = splitLine[0]
			splitInteraction = interaction.split("_") #first part is the enhancer, 2nd part the gene
			
			enhancerInfo = splitInteraction[0]
			splitEnhancerInfo = enhancerInfo.split(":")
			#Add the chr notation for uniformity. 		
			chrMatch = re.search("chr", splitEnhancerInfo[0], re.IGNORECASE)
			chrName = ""
			if chrMatch is None:
				chrName = "chr" + splitEnhancerInfo[0]
			else:
				chrName = splitEnhancerInfo[0]
			splitPosInfo = splitEnhancerInfo[1].split("-") #get the positions
			start = int(splitPosInfo[0])
			end = int(splitPosInfo[1])
			
			#Get the gene name
			splitGeneInfo = splitInteraction[1].split("$")
			geneName = splitGeneInfo[1]
			
			outF.write(chrName + "\t" + str(start) + "\t" + str(end) + "\t" + geneName + "\n")