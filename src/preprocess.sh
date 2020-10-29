### PREPROCESSING OF REGULATORY INFORMATION AND INPUT DATA###

# Some regulatory data files and other files need to be parsed to work properly as input to the tool.

#set 'false' to 'true' to run a specific block.

### PARSE CPG ISLANDS ###
run=false

if $run; then
	runFolder='./DataProcessing'
	inFile='../data/cpg/cpgIslandExt.txt'
	outFile='../data/cpg/cpg2.bed'

	python "$runFolder/cpgToBed.py" "$inFile" "$outFile"

fi

### CANCER TYPE-SPECIFIC DATA PROCESSING ###
#These are the cancer types we want to run the tool for, so preprocess their data.
types=('urinaryTract' 'prostate' 'esophagus' 'skin' 'pancreas' 'uterus' 'kidney' 'nervousSystem')


### PARSE EQTLS ###
#This requires download of GTEx_Analysis_v7.metasoft.txt.gz. See eQTLs/readme.txt.
run=false

if $run; then
	runFolder='./DataProcessing'
	inFile='../data/eQTLs/GTEx_Analysis_v7.metasoft.txt'
	geneLookupFile='../data/genes/ensemblGenesHg19'

	#loop over the cancer types
	for type in ${types[@]}; do
		outFolder='../data/eQTLs/$type/'
		python "$runFolder/filterEQTLs.py" "$inFile" "$geneLookupFile" "$outFolder" "$type"

		#for a few types, we do not have tissue-specific data, so we use blood.
		if [ "$type" = "urinaryTract" ] || [ "$type" = "kidney" ]; then

			python "$runFolder/filterEQTLs.py" "$inFile" "$geneLookupFile" "$outFolder" "blood"

		fi

	done

fi

run=false
if $run; then
	runFolder='./DataProcessing/'

	#loop over the different types
	for type in ${types[@]}; do

		#Cluster histones
		python "$runFolder/genericClustering.py" "../data/histones/$type/"

		#cluster eQTLs
		python "$runFolder/eQTLClustering.py" "../data/eQTLs/$type/"
	done

fi


### TMM NORMALIZATION ###
run=false

if $run; then
	runFolder='./DataProcessing'
	settingsFile='./settings/settings_HMF_Breast' #which file does not matter, all we need
	#are the paths to the svDir, expressionDir and metadata file.

	#first merge the HMF expression files
	python "$runFolder/mergeHMFExpression.py" "$settingsFile"
	
	Rscript "$runFolder/normalizeEdgeR.R"

fi

