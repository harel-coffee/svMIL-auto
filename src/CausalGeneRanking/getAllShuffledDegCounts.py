import sys
import numpy as np
import glob


def getAllCounts(files):
	
	#go through the files and get the number
	counts = []
	for currentFile in files:
		count = np.loadtxt(currentFile)
		counts.append(count)
	
	return counts	


shuffledPath = sys.argv[1]
# 
windowedCosmicCounts = getAllCounts(glob.glob(shuffledPath + 'windowedCosmicDegPairs.txt*'))
tadCosmicCounts = getAllCounts(glob.glob(shuffledPath + 'tadCosmicDegPairs.txt*'))
rulesCosmicCounts = getAllCounts(glob.glob(shuffledPath + 'rulesCosmicDegPairs.txt*'))

print("windowed cosmic+deg: ", np.mean(windowedCosmicCounts))
print("tad cosmic+deg: ", np.mean(tadCosmicCounts))
print("rules cosmic+deg: ", np.mean(rulesCosmicCounts))

# print "no of genes in the cosmic case for rules: ", len(ruleGenesCosmic)
# plt.hist(rulesCosmicCounts)
# plt.show()
# plt.clf()

windowedBcCounts = getAllCounts(glob.glob(shuffledPath + 'windowedBcDegPairs.txt*'))
tadBcCounts = getAllCounts(glob.glob(shuffledPath + 'tadBcDegPairs.txt*'))
rulesBcCounts = getAllCounts(glob.glob(shuffledPath + 'rulesBcDegPairs.txt*'))


print("windowed bc+deg: ", np.mean(windowedBcCounts))
print("tad bc+deg: ", np.mean(tadBcCounts))
print("rules bc+deg: ", np.mean(rulesBcCounts))

# print "no of genes in the cosmic case for rules: ", len(ruleGenesBc)
# plt.hist(rulesBcCounts)
# plt.show()
# plt.clf()

windowedDegCounts = getAllCounts(glob.glob(shuffledPath + 'windowedDegPairs.txt*'))
tadDegCounts = getAllCounts(glob.glob(shuffledPath + 'tadDegPairs.txt*'))
rulesDegCounts = getAllCounts(glob.glob(shuffledPath + 'rulesDegPairs.txt*'))

print("windowed deg: ", np.mean(windowedDegCounts))
print("tad deg: ", np.mean(tadDegCounts))
print("rules deg: ", np.mean(rulesDegCounts))

windowedCosmicCounts = getAllCounts(glob.glob(shuffledPath + 'windowedCosmicPairs.txt*'))
tadCosmicCounts = getAllCounts(glob.glob(shuffledPath + 'tadCosmicPairs.txt*'))
rulesCosmicCounts = getAllCounts(glob.glob(shuffledPath + 'rulesCosmicPairs.txt*'))

print("windowed cosmic: ", np.mean(windowedCosmicCounts))
print("tad cosmic: ", np.mean(tadCosmicCounts))
print("rules cosmic: ", np.mean(rulesCosmicCounts))

# print "no of genes in the cosmic case for rules: ", len(ruleGenesCosmic)
# plt.hist(rulesCosmicCounts)
# plt.show()
# plt.clf()

windowedBcCounts = getAllCounts(glob.glob(shuffledPath + 'windowedBcPairs.txt*'))
tadBcCounts = getAllCounts(glob.glob(shuffledPath + 'tadBcPairs.txt*'))
rulesBcCounts = getAllCounts(glob.glob(shuffledPath + 'rulesBcPairs.txt*'))


print("windowed bc: ", np.mean(windowedBcCounts))
print("tad bc: ", np.mean(tadBcCounts))
print("rules bc: ", np.mean(rulesBcCounts))

