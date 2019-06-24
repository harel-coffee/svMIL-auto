import sys
import re


with open(sys.argv[2], 'w') as outF:
	with open(sys.argv[1], 'r') as inF:
		
		for line in inF:
			if re.match("^#", line):
				outF.write(line)
				continue
			
			splitLine = line.split("\t")
			
			if splitLine[6] == "PASS":
				if splitLine[4] == "<DUP>" or splitLine[4] == "<DEL>" or splitLine[4] == "<CTX>" or splitLine[4] == "<INV>":
					outF.write(line)
		
		
		
		
		

