"""
	Some functions to plot the quality of the ranking. 

"""

#Plot the number of SVs and number of eQTLs by rank
#the X axis is the rank
#The number of SVs and number of eQTLs are on the y-axis

from __future__ import absolute_import
import sys
from six.moves import range

rankFile = sys.argv[1]

numberOfSVs = []
numberOfEQTLs = []

with open(rankFile, 'r') as inF:
	
	for line in inF:
		
		line = line.strip()
		splitLine = line.split("\t")
		
		#Get the number of SVs
		numberOfSVs.append(float(splitLine[2]))

		#Get the number of eQTLs
		numberOfEQTLs.append(float(splitLine[4]))
		

import matplotlib.pyplot as plt

x = list(range(0, len(numberOfSVs)))
#plt.scatter(x,numberOfSVs)

plt.scatter(x, numberOfEQTLs)

plt.show()
		
		
		







