from __future__ import absolute_import
import numpy as np
import pylab
import sys
from six.moves import range

dataDict = np.load(sys.argv[1]).item()
randomDataDict = np.load(sys.argv[2]).item()

data = []
for key in dataDict:
	for ind in range(0, dataDict[key]):
		data.append(key)

randomData = []
for key in randomDataDict:
	for ind in range(0, randomDataDict[key]):
		randomData.append(key)


#Calculate quantiles
data.sort()
quantile_levels1 = np.arange(len(data),dtype=float)/len(data)

randomData.sort()
quantile_levels2 = np.arange(len(randomData),dtype=float)/len(randomData)

#Use the smaller set of quantile levels to create the plot
quantile_levels = quantile_levels2

#We already have the set of quantiles for the smaller data set
quantiles2 = randomData

#We find the set of quantiles for the larger data set using linear interpolation
quantiles1 = np.interp(quantile_levels,quantile_levels1,data)

#Plot the quantiles to create the qq plot
pylab.scatter(quantiles1,quantiles2)
pylab.xlabel('True distribution')
pylab.ylabel('Random SV distribution')

#Add a reference line
maxval = max(data[-1],randomData[-1])
minval = min(data[0],randomData[0])
pylab.plot([minval,maxval],[minval,maxval],'k-')

pylab.show()