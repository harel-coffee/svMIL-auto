#!/usr/bin/env python

######### Documentation #########

__author__ = "Marleen Nieboer"
__credits__ = []
__maintainer__ = "Marleen Nieboer"
__email__ = "m.m.nieboer@umcutrecht.nl"
__status__ = "Development"

"""
	This is a very simple script that can read HiC interactions into a dictionary of measured interaction counts and their observed frequency in the file. These data can be plotted.
	With these plots, it helps to determine where to possibly set a threshold for the number of measured interactions before an interaction is considered confident. For example,
	interactions that have been counted only once may not be interesting or confident. The frequencies are log scaled in the plot since there are a lot of interactions that
	are measured only once.
	
	The first part of the code generates the dictionary, but it is a slow process, so it is recommended to run this once, let the dictionary be stored in a pickle file, which can
	then quickly be loaded back into memory by Python in the second part. The second part is used to do the plotting, which can be experimented with quickly as the files do not
	need to be re-read each time. 

"""

### Imports ###
import numpy as np
import matplotlib.pyplot as plt
import os
import re
import time
import glob
import pickle

#Uncomment this part of the code until the end of part 1 to read the measured interaction counts from the interaction matrices.



hiCFolder = "../../../../../Data/HiC/HUVEC_interchromosomal/5kb_resolution_interchromosomal/"

#Obtain the names of the files that we wish to read for each folder, use a wildcard since the names are different
filesToRead = []
for folderName in os.listdir(hiCFolder):
	#Skip if it is an osx specific file and not a folder
	if os.path.isdir(hiCFolder + "/" + folderName):
		filesToRead.append(hiCFolder + "/" + folderName + "/MAPQG0/*.RAWobserved")

startTime = time.time()

interactionCounts = dict()

fileCount = 0
for hiCFilePath in filesToRead:
	hiCFile = glob.glob(hiCFilePath) #get the actual file path instead of the wildcard one
	print hiCFile[0] #There will always be only one file that matches
	
	currentEdges = [] #refresh the list to make sure that not everything is in memory
	with open(hiCFile[0], "r") as f:
		for line in f:
			line = line.strip()
			splitLine = line.split("\t")
		
			count = int(float(splitLine[2]))
			if count > 0:
				
				if count in interactionCounts:
					interactionCounts[count] += 1
				else:
					interactionCounts[count] = 1
	

	fileCount += 1
	
endTime = time.time()
print "Took ", endTime - startTime, " seconds to store edges from all files"

#store the data with pickle, do not re-read the files each time

output = open('interactionCounts_interchromosomal.pkl', 'wb')

# Pickle dictionary using protocol 0.
pickle.dump(interactionCounts, output)

exit()


########## END PART 1




#Read the interactions back into memory
pkl_file = open('interactionCounts_intrachromosomal.pkl', 'rb')

interactionCounts = pickle.load(pkl_file)

#print some statistics of the interaction counts
print interactionCounts

print "highest interaction count: ", max(interactionCounts.keys())
print "highest interaction frequency: ", max(interactionCounts.values())

print "lowest interaction count: ", min(interactionCounts.keys())
print "lowest interaction frequency: ", min(interactionCounts.values())

#Here do the filtering to determine which values we want or not in our dictionary

filteredDictionary = dict()
threshold = 0
for interactionCount in interactionCounts:
	interactionValue = interactionCounts[interactionCount]
	
	if interactionCount > threshold:
		filteredDictionary[interactionCount] = np.log(interactionValue) #do a log to make the data more insightful
		
print filteredDictionary

#Do the plotting

plt.bar(filteredDictionary.keys(), filteredDictionary.values(), 50)
plt.xlabel("Measured interaction count")
plt.ylabel("Frequency (log)")
axes = plt.gca()
plt.show()
#plt.savefig('interchromosomal_allInteractions_lt10.svg')


