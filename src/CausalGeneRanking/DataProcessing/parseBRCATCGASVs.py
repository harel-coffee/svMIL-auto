

import sys

inFile = sys.argv[1]
outFile = sys.argv[2]

with open(outFile, 'w') as outF:
	outF.write("chr1\ts1\te1\to1\tchr2\ts2\te2\to2\tsource\tsample_name\tsv_type\tcancer_type\n")
	with open(inFile, 'r') as inF:
		
		lineCount = 0
		for line in inF:
			
			if lineCount < 2:
				lineCount += 1
				continue
			
			splitLine = line.split("\t")
			
			#Parse to: chr1	s1	e1	o1	chr2	s2	e2	o2	source	sample_name	sv_type	cancer_type
			sampleName = splitLine[0]
			chr1 = splitLine[4]
			s1 = splitLine[5]
			e1 = splitLine[5]
			o1 = splitLine[6]
			if o1 == "1":
				o1 = "+"
			if o1 == "-1":
				o1 = "-"
			
			chr2 = splitLine[9]
			s2 = splitLine[10]
			e2 = splitLine[10]
			o2 = splitLine[11]
			
			if o2 == "1":
				o2 = "+"
			if o2 == "-1":
				o2 = "-"
			
			svType = splitLine[14]
			
			outF.write(chr1 + "\t" + s1 + "\t" + e1 + "\t" + o1 + "\t" + chr2 + "\t" + s2 + "\t" + e2 + "\t" + o2 + "\t" + "Zhang et al 2018" + "\t" + sampleName + "\t" + svType + "\t" + "BRCA\n")
			
			
			

