"""
	Take SV file as input and filter for the given cancer type. 

"""

import re
import sys

inFile = sys.argv[1]
outFile = sys.argv[2]

cancerType = "ovarian"

with open(outFile, 'wb') as outF:
	with open(inFile, 'rb') as inF:
		lineCount = 0
		for line in inF:
			
			if lineCount < 1:
				outF.write(line)
				lineCount += 1
				continue
			
			line = line.strip()
			
			typeMatch = re.search(cancerType, line, re.IGNORECASE)
			if typeMatch is not None:
				outF.write(line + "\n")
			
		


