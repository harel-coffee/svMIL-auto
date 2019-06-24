"""
	Go through all the transcription factor files and make 1 big bed file

"""

import sys
import os
from os import listdir
from os.path import isfile, join

tfFolder = sys.argv[1]
outFile = sys.argv[2]

tfFiles = [f for f in listdir(tfFolder) if isfile(join(tfFolder, f))]

uniqueTFs = dict()

for tfFile in tfFiles:
	
	with open(tfFolder + tfFile, 'r') as inF:
		for line in inF:
			
			line = line.strip()
			splitLine = line.split("\t")
			
			chrName = splitLine[0]
			start = splitLine[1]
			end = splitLine[2]
			
			uuid = chrName + "\t" + start + "\t" + end
			if uuid not in uniqueTFs:
				uniqueTFs[uuid] = [chrName, start, end]
				
				
with open(outFile, 'w') as outF:
	
	for uniqueTF in uniqueTFs:
		
		outF.write(uniqueTF + "\n")
		
			
	