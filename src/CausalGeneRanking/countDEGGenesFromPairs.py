


from __future__ import absolute_import
from __future__ import print_function
import sys
import numpy as np

pairs = np.loadtxt(sys.argv[1], dtype="object")

genes = dict()
for pair in pairs[:,0]:
	
	splitPair = pair.split("_")
	genes[splitPair[0]] = 0
	
	
print(genes)
print(len(genes))