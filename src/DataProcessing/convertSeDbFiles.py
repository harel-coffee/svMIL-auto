"""
	Convert the bed files from SEDb to actual bed files

"""

import sys
import glob

inFolder = sys.argv[1]
inFile = glob.glob(inFolder + '*')[0]
outFile = sys.argv[2]

with open(outFile, 'w') as outF:
	
	with open(inFile, 'r') as inF:
	
		lineCount = 0
		for line in inF:
	
			line = line.strip()
			if lineCount < 1:
				lineCount += 1
				continue
	
			splitLine = line.split('\t')
	
			chrom = splitLine[2]
			start = splitLine[3]
			end = splitLine[4]
			
			outF.write(chrom + '\t' + start + '\t' + end + '\n')
			