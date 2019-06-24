"""
	Filter the COSMIC dataset to one specific type of cancer

"""

import sys
import re

cosmicFile = sys.argv[1]
cancerType = sys.argv[2]
outFile = sys.argv[3]

with open(outFile, 'w') as outF:
	with open(cosmicFile, 'rb') as inF:
		lineCount = 0
		for line in inF:
			
			if lineCount < 1:
				lineCount += 1
				continue
			
			
			
			splitLine = line.split("\t")
			typeMatch = re.search(cancerType, splitLine[9], re.IGNORECASE)
			
			if typeMatch is not None:
				outF.write(splitLine[0] + "\n")
			