#!/bin/bash
#$ -V
#$ -S /bin/bash
#$ -cwd
#$ -l h_vmem=16G
#$ -l h_rt=24:00:00
#$ -e workflow_err
#$ -o workflow_out

### INSTRUCTIONS ###

#	This script intends to show the workflow used to generate all figures in the paper.

#	The workflow is listed in order. 'Required' indicates a step that is needed to generate
#	output that is used in all figures. Only skip this steps if you have already generated
#	the data, and want to re-use that in a figure.

#	If you want to skip a part of the workflow, change the 'true' to 'false' for each step.

#	Memory requirements vary per workflow step and based on your dataset size.
#	For the HMF data, at least 16 GB of memory is required to load the bags for MIL in memory.

## workflow specific for the gnomAD test to see if the results are different from DGV

### (REQUIRED) PART 1 - DATA AND PATHS ###

#Load in a settings file with the paths to the data that all code will be run on.
## path here
settingsFolder='./settings/settings_HMF_BRCA/'

#Create a folder in which all output for this data will be stored
#Different steps will create their own intermediate folders in here
outputFolder='output/HMF_BRCA_gnomAD'

### (REQUIRED) PART 2 - LINK SVS TO GENES ###
run=false #Only skip this step if all output has already been generated!

if $run; then
	runFolder='./linkSVsGenes/'
	#Map the SVs to genes. This also outputs bags for MIL.
	python "$runFolder/main.py" "" "False" "0" "$settingsFolder" "$outputFolder"

fi

### FIGURE 2 - 2E: HEATMAP ###
run=true

if $run; then
	runFolder='./linkSVsGenes/'

	#Map the SVs to genes. Specific for germline, settings for this are different
	settingsFolder='./settings/settings_HMF_BRCA_gnomAD/'
	#python "$runFolder/main.py" "germline" "False" "0" "$settingsFolder" "$outputFolder"

	#Repeat for shuffled SVs
	#settingsFolder='./settings/settings_HMF_BRCA/'
	#python "$runFolder/main.py" "random" "True" "0" "$settingsFolder" "$outputFolder"

	#split affected/non-affected pairs
	#python "$runFolder/splitPairsPathogenicNonPathogenic.py" "$outputFolder"

	#Then make the heatmap plot
	python "$runFolder/plotPathogenicNonPathogenicFeatures.py" "$outputFolder"

fi
