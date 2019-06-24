import sys
import numpy as np
import os
from os import listdir
from os.path import isfile, join

snvFolder = sys.argv[1]

snvFiles = [f for f in listdir(snvFolder) if isfile(join(snvFolder, f))]

genes = []
for snvFile in snvFiles:
	if snvFile == "MANIFEST.txt":
		continue
	with open(snvFolder + snvFile, 'r') as inF:
		lineCount = 0
		for line in inF:
			if lineCount < 1:
				lineCount += 1
				continue
			
			line = line.strip()
			splitLine = line.split("\t")
			
			if splitLine[25] == "Somatic":
				genes.append(splitLine[0])

print np.unique(genes)
print len(np.unique(genes))
