# svMIL2: predicting the pathogenic effect of somatic non-coding structural variants disrupting the 3D genome through Multiple Instance Learning

svMIL2 is the improved version of svMIL, tested across multiple cancer types.

For details on the previous version of svMIL, refer to: https://doi.org/10.1093/bioinformatics/btaa802

# How to use

An example test dataset is provided in the src/test folder.

svMIL requires the following data for each cancer patient: SVs, SNVs, CNVs and expression data. Patients
without CNVs, SNVs or expression data will be skipped automatically. These data are expected in the
following folder structure to make it compatible with the HMF data used in the paper:

somatics/ (here CNVs, SNVs and SV data should be stored)
	- somatics_patient1 (make sure that there is an underscore before the patient ID)
		- patient1.sv.ann.vcf.gz (SV data in VCF format, this exact name is required)
		- patient1.somatic.vcf.gz (SNV calls in VCF format, this exact name is required)
		- patient1.cnv.gene.tsv (CNV calls in tsv format, listing the copy numbers of each gene, this exact name is required)
	- somatics_patient2
		- ...
expression/ (here, expression data (read counts) per patient should be stored)
	- patient1 (make sure that the patient IDs match the patient IDs in the somatics folder, the folder names are used to check if expression data exists for each patient)
		- read counts file (not used by the model directly, see expression data section below)
	- patient2


## Step 1: setting up the regulatory data (features)

Used regulatory data is provided in the data folder. For more details on data sources and used tissue types, please refer to Table S1 (TBA).

Most required pre-processed data is provided in the data folder. Only 2 things are needed:

- eQTLs are too large to provide and need to be parsed and clustered.
- Normalization of expression data (see expression data section below)

## Step 2: formatting input data

Examples of the required files are in the src/test folder.

#SV calls

SV calls are expected in gzipped VCF format. The CHROM, POS, ALT and INFO fields are required. The INFO field
is used to only filter out calls with a PASS. From the ALT field, the second chromosome, position and SV
orientation are obtained.

#SNV calls

SNV calls are expected in gzipped VCF format. The CHROM, POS, INFO and FORMAT fields are required. The INFO field
is used to filter out calls with a PASS. From the FORMAT field, snpEff annotated gene names are parsed
to link SNVs to genes (which should be formatted like TF_binding_site_variant|MODIFIER|GJB3)

#CNV calls

CNV calls are expected in tsv format. The chromosome, start, end, gene, minCopyNumber and maxCopyNumber fields
are required. The copy number of each gene is determined through this file. Genes with a copy number
< 1.7 or > 2.3 are considered affected by a copy number. 

#Expression data

svMIL2 requires a normalized expression data file with the patient IDs in the columns, and gene
names in the rows. For the HMF data, we used TMM normalization.

#Metadata

The metadata file is used to link patient IDs to cancer type. This is useful if the data folder
contains data of multiple cancer types, and the model should be run on only 1 type. The columns
patientId, sampleId and primaryTumorLocation (cancer type) are required. Patient ID and sample ID
should be the same.

## Step 3: setting paths

All paths to data and data folders should be in a file called settings.py.

## Step 4: running svMIL2

All commands required to run the model are provided in src/test/test.sh. This script should be
run from the main src folder and provide the path to the settings file folder and the output folder, e.g.:

sh test/test.sh test/settings test/output

Please refer to test.sh for more details.

## Step 5: output files

Within the output folder, a number of files are created that may be relevant to the user.

linkedSVGenePairs
	- nonCoding_geneSVPairs.txt (List of all SV-gene pairs and which regulatory elements/features are disrupted)
	- nonCoding_geneSVPairs.txt_pathogenicPairsFeatures.txt & nonCoding_geneSVPairs.txt_nonPathogenicPairsFeatures.txt (SV-gene pairs split into pathogenic and non-pathogenic based on z-score)
	- bags.pkl & normalizedBags.pkl (bags made by svMIL)
tadDisruptionsZScores
	- zScores.txt (For each gene in each patient, the z-score from patients with TAD disruptions to patients without TAD disruptions. Only contains pairs WITHOUT coding mutations)
multipleInstanceLearning
	- leaveOnePatientOutCV
		- leaveOnePatientOutCV_SVType.txt (SV-gene pair, true label, predicted label)
rocCurves
	- Output ROC and PR curves from the leave-one-patient-out CV run.



# How to use: recreating paper figures

## Step 1: setting up the data

All used data in the paper (SVs, SNVs, CNVs and RNA-seq data) needs to be requested from the HMF: https://www.hartwigmedicalfoundation.nl/en/appyling-for-data/.

## Step 2: preprocessing of the data

The TMM normalization step of the RNA-seq data from the HMF dataset is required. Other data needs
to be pre-processed as detailed in step 1 above.

## Step 3: generating all paper figures

plotting.sh lists all steps and instructions to generate the figures.


# Requirements

All code was tested using:
- Python 3.7.3
- Numpy 1.16.4
- Scipy 1.3.0
- Scikit-learn 0.22.2.postl
- Matplotlib 3.2.1
- Statsmodels 0.10.0


