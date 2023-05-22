#!/bin/bash

#################
#set a job name
#SBATCH --job-name=WGSassign
#SBATCH --output=wgsa.out
#SBATCH --error=wgsa.err
#################
#SBATCH --time=5:00:00
#SBATCH --ntasks-per-node=24
#################

source ~/.bashrc

conda activate WGSassign

BEAGLE="/WGSA_data/yewa.known.ind105.ds_2x.beagle.gz"
POP_ID="/WGSA_data/yewa.known.ind105.reference.IDs.txt"

WGSassign --beagle $BEAGLE --pop_af_IDs $POP_ID --get_reference_af --loo --out yewa.wgs.assignment --threads 20
