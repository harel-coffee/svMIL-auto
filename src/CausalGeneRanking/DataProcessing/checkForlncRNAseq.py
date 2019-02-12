
import sys

inFile = sys.argv[1]

lncChr = "chr17"
lncPos1 = 41447213
lncPos2 = 41466266

with open(inFile, 'r') as inF:
	lineCount = 0
	for line in inF:
		if lineCount < 2:
			lineCount += 1
			continue
		
		splitLine = line.split("\t")
		exon = splitLine[0]
		
		splitExonChr = exon.split(":")
		chrom = splitExonChr[0]
		
		splitExonPos = splitExonChr[1].split("-")
		pos1 = int(splitExonPos[0])
		pos2 = int(splitExonPos[1])
		
		if pos1 > pos2:
			tmpPos2 = pos2
			pos2 = pos1
			pos1 = tmpPos2
		
		# print chrom
		# print pos1
		# print pos2
		#
		
		#print "exon len: ", pos2 - pos1
		# print "lnc len: ", lncPos2 - lncPos1
		# exit()
		
		if chrom == lncChr:
			if lncPos1 < pos1 and lncPos2 > pos2:
				print "exon: ", splitLine[0]
		
		
				