
def getAllCounts(files):
	
	#go through the files and get the number
	counts = []
	for currentFile in files:
		count = np.loadtxt(currentFile)
		counts.append(count)
	
	return counts	


import glob

shuffledPath = sys.argv[1]
# 
windowedCosmicCounts = getAllCounts(glob.glob(shuffledPath + 'windowedCosmicDeg.txt*'))
tadCosmicCounts = getAllCounts(glob.glob(shuffledPath + 'tadCosmicDeg.txt*'))
rulesCosmicCounts = getAllCounts(glob.glob(shuffledPath + 'rulesCosmicDeg.txt*'))

print("windowed cosmic+deg: ", np.mean(windowedCosmicCounts))
print("tad cosmic+deg: ", np.mean(tadCosmicCounts))
print("rules cosmic+deg: ", np.mean(rulesCosmicCounts))

# print "no of genes in the cosmic case for rules: ", len(ruleGenesCosmic)
# plt.hist(rulesCosmicCounts)
# plt.show()
# plt.clf()

windowedBcCounts = getAllCounts(glob.glob(shuffledPath + 'windowedBcDeg.txt*'))
tadBcCounts = getAllCounts(glob.glob(shuffledPath + 'tadBcDeg.txt*'))
rulesBcCounts = getAllCounts(glob.glob(shuffledPath + 'rulesBcDeg.txt*'))


print("windowed bc+deg: ", np.mean(windowedBcCounts))
print("tad bc+deg: ", np.mean(tadBcCounts))
print("rules bc+deg: ", np.mean(rulesBcCounts))

# print "no of genes in the cosmic case for rules: ", len(ruleGenesBc)
# plt.hist(rulesBcCounts)
# plt.show()
# plt.clf()

windowedDegCounts = getAllCounts(glob.glob(shuffledPath + 'windowedDeg.txt*'))
tadDegCounts = getAllCounts(glob.glob(shuffledPath + 'tadDeg.txt*'))
rulesDegCounts = getAllCounts(glob.glob(shuffledPath + 'rulesDeg.txt*'))

print("windowed deg: ", np.mean(windowedDegCounts))
print("tad deg: ", np.mean(tadDegCounts))
print("rules deg: ", np.mean(rulesDegCounts))

