import sys
import re

cpgFile = sys.argv[1]
outFile = sys.argv[2]

with open(outFile, 'w') as outF:
	with open(cpgFile, 'r') as f:
			
		lineCount = 0
		for line in f:
			if lineCount < 1:
				lineCount += 1
				continue
			
			line = line.strip()
			splitLine = line.split("\t")
			
			#Add the chr notation for uniformity. 		
			chrMatch = re.search("chr", splitLine[1], re.IGNORECASE)
			chrName = ""
			if chrMatch is None:
				chrName = "chr" + splitLine[1]
			else:
				chrName = splitLine[1]
				
			start = int(splitLine[2])
			end = int(splitLine[3])
			
			outF.write(chrName + "\t" + str(start) + "\t" + str(end) + "\n")