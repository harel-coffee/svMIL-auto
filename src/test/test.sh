### INSTRUCTIONS ###

#	This workflow is intended to do a test run for svMIL2 with test data. It lists the steps
#	that would be part of a typical run with svMIL2.

#	If you want to skip a part of the workflow (e.g. to re-do one step only),
#	change the 'false' to 'true' for each step.

### (REQUIRED) PART 1 - DATA AND PATHS ###
#Load in a settings file with the paths to the data that all code will be run on.
## The path should be to the folder in which the settings.py file will be loaded e.g.
## test/settings
settingsFolder="$1"
#Create a folder in which all output for this data will be stored
#Different steps will create their own intermediate folders in here
## This is e.g. test/output
outputFolder="$2"

### (REQUIRED) PART 2 - LINK SVS TO GENES ###
run=false #Only skip this step if all output has already been generated!

if $run; then
	runFolder='./linkSVsGenes/'
	#Map the SVs to genes. This also outputs bags for MIL.
	python "$runFolder/main.py" "" "False" "0" "$settingsFolder" "$outputFolder"

fi

### (REQUIRED) PART 3 - IDENTIFY PATHOGENIC SV-GENE PAIRS ###
run=false

if $run; then
	runFolder='./tadDisruptionsZScores/'

	#first link mutations to patients. These are required to quickly check which patients
	#have which mutations in which genes
	python "$runFolder/determinePatientGeneMutationPairs.py" "$settingsFolder" "$outputFolder"

	#identify which TADs are disrupted in these patients, and compute the
	#z-scores of the genes in these TADs. Filter out genes with coding mutations.
	python "$runFolder/computeZScoresDisruptedTads.py" "$settingsFolder" "$outputFolder" "False"

	#split the SV-gene pairs into pathogenic/non-pathogenic, which we use later on.
	runFolder='./linkSVsGenes/'
	python "$runFolder/splitPairsPathogenicNonPathogenic.py" "$outputFolder"

fi

### LEAVE-ONE-PATIENT-OUT CV SCENARIO ###

### PART 1 - SETTING UP FOR MULTIPLE INSTANCE LEARNING ###
run=false

if $run; then
	runFolder='./multipleInstanceLearning/'

	#first normalize the bags
	python "$runFolder/normalizeBags.py" "$outputFolder"

	#then generate the similarity matrices for all SVs in the lopoCV case
	python "$runFolder/generateSimilarityMatrices.py" "$outputFolder" "False" "True" "False" "False" "False" "$settingsFolder"

fi

### PART 2 - LEAVE-ONE-PATIENT-OUT CV FOR MIL PERFORMANCE ###
run=false

if $run; then
	runFolder='./multipleInstanceLearning/'

	#test the classifier and output the MIL curves
	python "$runFolder/runMILClassifier.py" "$outputFolder" "False" "True" "False" "False" "False"

fi


### TRAIN ON ONE DATASET, APPLY TO ANOTHER SCENARIO ###

#For this example, we have only one dataset. So we train/test on the same data to illustrate.
#Normally, if you have 2 datasets, you would run step 1-3 above for both of those datasets, and
#provide the output folder with the training dataset here in step 1.

trainOutputFolder='test/output'
testOutputFolder='test/output'

### PART 1 - GENERATE SIMILARITY MATRIX ON TRAINING DATA ###
run=true

if $run; then
	runFolder='./multipleInstanceLearning/'

	#first normalize the bags
	python "$runFolder/normalizeBags.py" "$trainOutputFolder"

	#This command generates a FULL similarity matrix on all patients, so not a leave-one-patient-out CV
	python "$runFolder/generateSimilarityMatrices.py" "$trainOutputFolder" "False" "False" "False" "False" "True" "$settingsFolder"

fi

#Then for step 2, we provide the output locations of both the training and test data to train
#the classifier. This script will also output a prioritized list of SV-gene pairs by pathogencity.

### PART 2 - ###
run=true

if $run; then
	runFolder='./multipleInstanceLearning/'

	#also normalize the bags for the test data.
	python "$runFolder/normalizeBags.py" "$testOutputFolder"

	python "$runFolder/applyClassifier.py" "$trainOutputFolder" "$testOutputFolder"
	#python "$runFolder/runMILClassifier.py" "$outputFolder" "False" "True" "False" "False" "False"

fi
