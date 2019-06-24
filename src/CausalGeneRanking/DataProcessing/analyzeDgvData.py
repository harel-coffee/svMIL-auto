import sys
import numpy as np
import re
# 
# dgvVariants = np.loadtxt(sys.argv[1], dtype="object")
# dgvVariants = dgvVariants[1:,:]

svTypes = []

with open(sys.argv[1], 'r') as inF:
	lineCount = 0
	for line in inF:
		if lineCount < 1:
			lineCount += 1
			continue
		line = line.strip()
		splitLine = line.split("\t")
		svType = splitLine[10]
		if svType != "del" and svType != "invers" and svType != "tandem_dup":
					
			interChrTypeMatch = re.search("chr", svType, re.IGNORECASE)
			transTypeMatch = re.search("trans", svType, re.IGNORECASE)
			rangeTypeMatch = re.search("range", svType, re.IGNORECASE)
			if interChrTypeMatch is None and transTypeMatch is None and rangeTypeMatch is None:
				continue

		
		svTypes.append(splitLine[10])
	
	


import matplotlib.pyplot as plt

#plt.hist(dgvVariants[:,10])
plt.hist(svTypes)
plt.show()


