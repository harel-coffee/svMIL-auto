import sys

geneScoreFile = sys.argv[1]
hotNetOutFile = sys.argv[2]
hotNetOutFile2 = sys.argv[3]


with open(hotNetOutFile, 'w') as outF1:
	with open(hotNetOutFile2, 'w') as outF2:
		with open(geneScoreFile, 'r') as inF:
			
			for line in inF:
				line = line.strip()
				splitLine = line.split("\t")
				
				outF1.write(splitLine[0] + "\t" + splitLine[1] + "\n")
				outF2.write(splitLine[0] + "\t" + splitLine[4] + "\n")
			
		