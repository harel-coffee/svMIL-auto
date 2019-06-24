"""
	Script to parse the oreganno TF binding sites. 

"""

import sys
import numpy as np

inFile = sys.argv[1]
outFile = sys.argv[2]

tfs = []
with open(outFile, 'w') as outF:
	with open(inFile, 'r') as inF:
		lineCount = 0 
		for line in inF:
			if lineCount < 1:
				lineCount += 1
				continue
			
			
			line = line.strip()
			splitLine = line.split("\t")
			
			#get the  chr, start and end.
			#filter for hg19 (13) and by type (3)
			
			
			if splitLine[3] == "TRANSCRIPTION FACTOR BINDING SITE" and splitLine[13] == "hg19":
				tfs.append([splitLine[15], splitLine[16], splitLine[17], splitLine[3], splitLine[13]])
				newLine = splitLine[15] + "\t" + splitLine[16] + "\t" + splitLine[17] + "\n"
				outF.write(newLine)
				
tfs = np.array(tfs)
print tfs
print tfs.shape
