import sys

indexFile = sys.argv[1]
outFile = sys.argv[2]

with open(outFile, 'w') as outF:
	with open(indexFile, 'r') as inF:
		
		for line in inF:
			
			
			splitLine = line.split("\t")
			splitSplitLine = splitLine[0].split(" ")
			outF.write(splitSplitLine[0] + "\t" + splitSplitLine[1] + "\n")