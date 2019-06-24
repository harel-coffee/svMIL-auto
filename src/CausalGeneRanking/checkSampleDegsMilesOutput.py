

import sys
import numpy as np

milesGenes = np.loadtxt(sys.argv[1], dtype="object")
degGenes = np.loadtxt(sys.argv[2], dtype="object")


print "Intersect between MILES genes and DEGs: ", np.intersect1d(milesGenes, degGenes[:,0])
print "Difference between MILES genes and DEGs: ", np.setdiff1d(milesGenes, degGenes[:,0])

print "Intersect size between MILES genes and DEGs: ", len(np.intersect1d(milesGenes, degGenes[:,0]))
print "Difference size between MILES genes and DEGs: ", len(np.setdiff1d(milesGenes, degGenes[:,0]))


