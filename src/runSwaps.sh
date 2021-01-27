#Run svMIL2 with swaps, and wait each time until another job has finished.
#To run for only one type, replace types with an array with only 1 type.

#We use this to call the right workflow.
#So e.g. if we want to start from breast and do all runs with data from other
#cancer types, ownType should be 'Breast', matching the file names in workflows.
ownType="$1"

#1. Generate the settings files automatically
python createSettingsFiles.py "$ownType"

#2. Do the runs
types=('hmec' 'ov' 'gm12878' 'coad' 'luad' 'urinaryTract' 'prostate' 'esophagus' 'skin' 'pancreas' 'uterus' 'nervousSystem' 'kidney')


RES='None'
for type in ${types[@]}; do

	if [ "$type" = "$ownType" ]; then
		continue
	fi

	runFolder='./workflows/'
	settingsFile="./settings/settings_HMF_${ownType}_${type}"
	outputFile="./output/HMF_${ownType}_${type}"

	#perform run with first type
	if [ "$RES" = "None" ]; then
		RES=$(sbatch --parsable "$runFolder/"workflow_"$ownType".sh "$settingsFile" "$outputFile")

		continue
	fi

	#only do the next run if the previous one has finished, to avoid creating
	#a lot of tmp space at once
	RES=$(sbatch -d afterany:${RES} --parsable "$runFolder/"workflow_"$ownType".sh "$settingsFile" "$outputFile")

done
