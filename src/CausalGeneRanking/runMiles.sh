#!/bin/bash
#$ -V
#$ -S /bin/bash
#$ -cwd
#$ -m as
#$ -M m.m.nieboer@umcutrecht.nl
#$ -l h_vmem=50G
#$ -l h_rt=48:00:00
#$ -e miles_lassoAcc2patients_err
#$ -o miles_lassoAcc2patients_out

python milesAllBags.py