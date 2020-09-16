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

### PARSE EQTLS ###
#This requires download of GTEx_Analysis_v7.metasoft.txt.gz. See eQTLs/readme.txt.
run=true

if $run; then
	runFolder='./DataProcessing'
	inFile='../data/eQTLs/GTEx_Analysis_v7.metasoft.txt'
	geneLookupFile='../data/genes/ensemblGenesHg19'

	breastOutFolder='../data/eQTLs/breast/'
	python "$runFolder/filterEQTLs.py" "$inFile" "$geneLookupFile" "$breastOutFolder" "breast"

	#repeat for ovary and liver
	ovarianOutFolder='../data/eQTLs/ov/'
	python "$runFolder/filterEQTLs.py" "$inFile" "$geneLookupFile" "$ovarianOutFolder" "ovarian"

	liverOutFolder='../data/eQTLs/liver/'
	python "$runFolder/filterEQTLs.py" "$inFile" "$geneLookupFile" "$liverOutFolder" "liver"

	liverOutFolder='../data/eQTLs/gm12878/'
	python "$runFolder/filterEQTLs.py" "$inFile" "$geneLookupFile" "$liverOutFolder" "blood"


fi

### PARSE HIC DATA ###
run=false

#Requires download of 5kb resolution intrachromosomal interactions. See hic/readme.txt

if $run; then
	runFolder='./DataProcessing'
	inFiles='../data/hic/HMEC/5kb_resolution_intrachromosomal/'
	thresholdedDataOutFolder='../data/hic/HMEC/thresholdedHiC/'

	#First apply a threshold to filter out regions with too few interactions.
	python "$runFolder/thresholdHiCData.py" "$inFiles" "$thresholdedDataOutFolder"

	#Then parse the thresholded interaction to 1 file
	parsedOutFile='../data/hic/HMEC_interactions.txt'
	python "$runFolder/parseHicInteractions.py" "$thresholdedDataOutFolder" "$parsedOutFile"

	#make tmp file to check if the same
	#Then group the interactions by the TAD that these take place in for faster processing
	groupedInteractionsFile='../data/hic/HMEC_groupedTADInteractions.txt'
	tadFile='../data/tads/hmec/HMEC_Lieberman-raw_TADs.bed'
	python "$runFolder/makeTADHiCFile.py" "$parsedOutFile" "$tadFile" "$groupedInteractionsFile"

fi

### PARSE GERMLINE VARIANTS ###
#Parse gnomAD variants. See svs/readme.txt
run=false
if $run; then
	runFolder='./DataProcessing'
	inFile='../data/svs/gnomad_v2.1_sv.sites.bed'
	outFile='../data/svs/gnomad_v2.1_sv.sites_filtered_01072020.bed'

	python "$runFolder/filterGnomadVariants.py" "$inFile" "$outFile"

fi

#Alternative to use DGV variants.
run=false
#Requires download of germline DGV variants. See svs/readme.txt

if $run; then
	runFolder='./DataProcessing'
	inFile='../data/svs/GRCh37_hg19_variants_2020-02-25.txt'
	outFile='../data/svs/GRCh37_hg19_variants_2020-02-25_filtered.txt'

	python "$runFolder/filterDgvVariants.py" "$inFile" "$outFile"

fi


### TMM NORMALIZATION ###
run=false

if $run; then
	runFolder='./DataProcessing'
	
	Rscript "$runFolder/normalizeEdgeR.R"

fi

### PROCESS GTEX DATA ###


### BATCH CORRECTION AND TMM NORMALIZATION OF BRCA AND GTEX EXPRESSION ###
run=false

if $run; then
	runFolder='./DataProcessing'
	gtexExpression='../data/expression/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_reads.gct'
	gtexMetaData='../data/expression/GTEx_v7_Annotations_SampleAttributesDS.txt'
	brcaExpressionFolder='/hpc/compgen/users/mnieboer/data/gtexPipeline/gtex/gtex_output/read_counts/' #this needs to be made
	
	#merge the BRCA read counts into one file, same format as the GTEx expression
	combinedReadCountsBRCAFile='../data/expression/combinedReadCountsBRCA.txt'
	#python "$runFolder/combineReadCountsBRCA.py" "$brcaExpressionFolder" "$combinedReadCountsBRCAFile"
	
	#### !needs proper output folder provided to write to
	#python "$runFolder/mergeBRCAGTEx.py" "$gtexExpression" "$gtexMetaData" "$combinedReadCountsBRCAFile" 
	
	#Then perform batch correction and TMM normalization
	### !also provide locations and output folders
	module load R/3.2.2
	Rscript "$runFolder/normalizeBRCAGTEx.R"
	
fi



### ANNOTATE HMF SVS WITH SV TYPES ###
run=false

if $run; then
	runFolder='./DataProcessing'

	Rscript "$runFolder/simple-event-annotation.R"

fi


### pre-process CNVs of TCGA ###
run=false

if $run; then
	runFolder='./DataProcessing'

	### TO DO:
	# -define out file
	# - make this run on all files, not specific to 1 cancer type
	# - maybe move it into the other pipeline to make it specific for the cancer type

	python "$runFolder/preprocessCNVs.py" "$settingsFolder" "$outFile"
fi
