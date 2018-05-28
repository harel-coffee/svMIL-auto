# Causal gene-based ranking

This branch is intended for the code where we rank known causal genes by how likely these were causal for a specific set of SVs and/or SNVs from cancer patients. An idea of the causal SVs and SNVs can then be obtained from the highest ranked genes. 

The currently included datasets are TADs and eQTLs. 

In the settings file the required file paths can be specified. 

# How to use

The starting script is runRankingWithPermutations.sh. This script does not require any parameters. If run on the HPC, it will first score all causal genes for causality and then repeat the process 1000 times with SVs and/or SNVs permuted across the genome. 

If these scripts are done, the computePValuesPerGene.py script can be run to do the actual ranking of the genes. As parameters the script requires the output folder containing the scores for the normal run and 1000 permutations, and the number of permutations that were run (+1 because there is currently still a bug :)).

The output is a list of all causal genes that have significant p-values in as much layers as possible. 

