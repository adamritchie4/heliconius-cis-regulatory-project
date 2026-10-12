#!/bin/bash

#SBATCH -J grelu
#SBATCH --ACCOUNT

#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --mem=10G
#SBATCH --time=3-00:00:00
#SBATCH --output=./logs/02.slurm-%j.out
#SBATCH --error=./logs/02.slurm-%j.err

##########################################################################################################
########################### Label setup and script to train gReLU CNN models #############################
##########################################################################################################

# Inputs:
#   - Training, test, and validation BED files
#
# Outputs:
#   - Merged labeled enhancer and background sequence TSV inputs
#   - Best gReLU model hyperparameters
#   - Trained gReLU dilated CNN

set -euo pipefail

##########################################################################################################
# Environment
##########################################################################################################

CONDA_BIN=$(dirname $(which conda))
CONDA_BASE=$(dirname $CONDA_BIN)

export PATH=$CONDA_BIN:$PATH

source $CONDA_BASE/etc/profile.d/conda.sh
conda activate heliconius_cnn_gpu

PYTHON=$(which python)

module load bedtools2/2.31.1-extr
module load samtools

mkdir -p $HOME/.config/genomepy

echo "CONDA_EXE=$CONDA_EXE"
echo "CONDA_PREFIX=$CONDA_PREFIX"
echo "Python=$PYTHON"
echo "Node=$(hostname)"
echo "SLURM_JOB_ID=${SLURM_JOB_ID:-unset}"
echo "CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES:-unset}"


##########################################################################################################
# gReLU setup
##########################################################################################################

PROJECT_ROOT=$(pwd)
CONFIG=$PROJECT_ROOT/config/config.yaml

echo "Project root: $PROJECT_ROOT"
echo "Config:       $CONFIG"

##########################################################################################################
# Training datasets
##########################################################################################################

HHIM_BEDS=$PROJECT_ROOT/data/training_species/Hhim/grelu
HLAT_BEDS=$PROJECT_ROOT/data/training_species/Hlat/grelu

HHIM_FASTA=$PROJECT_ROOT/$(yq -r '.training_species.Hhim.genome.fasta' $CONFIG)
HLAT_FASTA=$PROJECT_ROOT/$(yq -r '.training_species.Hlat.genome.fasta' $CONFIG)

TRAINING=$PROJECT_ROOT/data/training
mkdir -p $TRAINING

echo "Hhim BEDs:  $HHIM_BEDS"
echo "Hlat BEDs:  $HLAT_BEDS"
echo "Hhim FASTA: $HHIM_FASTA"
echo "Hlat FASTA: $HLAT_FASTA"
echo "Output:     $TRAINING"

##########################################################################################################

$PYTHON "$PROJECT_ROOT/src/prepare_grelu_data.py" \
   	--hhim_beds "$HHIM_BEDS" \
    	--hlat_beds "$HLAT_BEDS" \
    	--hhim_fasta "$HHIM_FASTA" \
   	--hlat_fasta "$HLAT_FASTA" \
    	--output "$TRAINING"

$PYTHON $PROJECT_ROOT/src/tune_grelu.py \
    	--train "$TRAINING/train.tsv" \
    	--valid "$TRAINING/valid.tsv" \
    	--output "$PROJECT_ROOT/models/grelu_tuning" \
    	--trials 20

$PYTHON $PROJECT_ROOT/src/train_grelu.py \
    	--train "$TRAINING/train.tsv" \
    	--valid "$TRAINING/valid.tsv" \
    	--test "$TRAINING/test.tsv" \
    	--output "$PROJECT_ROOT/models/grelu"
