# svMIL: predicting the pathogenic effect of somatic structural variants through multiple instance learning

# How to use: setting up the data

For the HMF data, SVs, SNVs, CNVs and RNA-seq data needs to be requested.

To run for the PCAWG data, make sure all datasets are properly downloaded. Some files are provided, but some will need to be downloaded as they are too large. The data folders that are incomplete and need to be downloaded are:

- data/cnvs
- data/expression
- data/snvs
- data/svs

Instructions are found in the readme.txt in the respective folders.

For the regulatory data, everything is provided except for the HiC data, as these files are too large. Follow the steps in data/hic/readme.txt to download this data.

All sources of the other files are also listed in the readme.txt files if these need to be re-downloaded.

# How to use: preprocessing

Except for Hi-C data, the provided data are already pre-processed. In the script src/preprocess.sh, steps are listed that are required to process the Hi-C data. Here, the processing steps for all
the other data is also listed, if needed.

# How to use: generating all paper figures

The script src/workflow.sh lists all steps that need to be executed to produces all figures in order. These steps are specific for the HMF breast cancer data.
The scripts src/workflow_pcawg_liver.sh and src/workflow_pcawg_ov.sh are specific for running on the PCAWG ovarian and liver datasets.

# Requirements

All code was tested using:
- Python 3.7.3
- Numpy 1.16.4
- Scipy 1.3.0
- Scikit-learn 0.22.2.postl
- Matplotlib 3.1.0
- Statsmodels 0.10.0


