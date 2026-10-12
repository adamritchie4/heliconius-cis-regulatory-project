#!/bin/bash

#SBATCH -J 01_prepare_labels
#SBATCH --account=BISC020662

#SBATCH --nodes=1
#SBATCH --cpus-per-task=5
#SBATCH --mem=10G
#SBATCH --time=1-00:00:00
#SBATCH --output=./logs/01.slurm-%A_%a.out
#SBATCH --error=./logs/01.slurm-%A_%a.err
#SBATCH --array=0-1

##########################################################################################################
#################### Prepare CNN labels and split into train/validation/test datasets ####################
##########################################################################################################

# Inputs:
#   - Enhancer BED files defined using H3K27ac and ATAC peaks
#   - TSS annotation BED files
#   - Genome FASTA assemblies
#   - Configurable paths from config/config.yaml
#
# Outputs:
#   - Processed enhancer BED files
#   - N-matched background BED files
#   - Training, validation, and test BED files for gReLU

set -euo pipefail

# input files are configurable from Heliconius_CisConservation_CNN/config/config.yaml and should be placed in HeliconiusCisConservation_CNN/data/training_species

##########################################################################################################
# Environment
##########################################################################################################

CONDA_BIN=$(dirname $(which conda))
export PATH=$CONDA_BIN:$PATH

CONDA_BASE=$(dirname $CONDA_BIN)
source $CONDA_BASE/etc/profile.d/conda.sh

conda activate heliconius_cnn

mkdir -p $HOME/.config/genomepy

echo "CONDA_EXE=$CONDA_EXE"
echo "CONDA_PREFIX=$CONDA_PREFIX"
echo "PATH=$PATH"

PYTHON=$(which python)

echo "Using Python: $PYTHON"

module load bedtools2/2.31.1-extr
module load samtools

##########################################################################################################
# Project setup
##########################################################################################################

echo "Preparing CNN training labels and gReLU datasets..."

PROJECT_ROOT=$(pwd)
CONFIG=$PROJECT_ROOT/config/config.yaml
SPLITS=$PROJECT_ROOT/config/splits.yaml

##########################################################################################################
# Select species from SLURM array
##########################################################################################################

mapfile -t SPECIES < <(yq -r '.training_species | keys | .[]' $CONFIG)

if [[ -z ${SLURM_ARRAY_TASK_ID:-} ]]; then
    echo "SLURM_ARRAY_TASK_ID not set. Running locally for species index 0."
    task_id=0
else
    task_id=$SLURM_ARRAY_TASK_ID
fi

if (( task_id >= ${#SPECIES[@]} )); then
    echo "ERROR: task_id=$task_id but only ${#SPECIES[@]} species are defined."
    exit 1
fi

sp=${SPECIES[$task_id]}

echo "Running pipeline for $sp"

##########################################################################################################
# Input paths
##########################################################################################################

GENOME=$PROJECT_ROOT/$(yq -r ".training_species.$sp.genome.fasta" $CONFIG)
ANNOTATION=$PROJECT_ROOT/$(yq -r ".training_species.$sp.annotation.tss" $CONFIG)
ENHANCERS=$PROJECT_ROOT/$(yq -r ".training_species.$sp.labels.enhancers" $CONFIG)

##########################################################################################################
# Output paths
##########################################################################################################

BASE=$PROJECT_ROOT/data/training_species/$sp

GENOME_O=$BASE/genome
LABELS=$BASE/labels
GRELU=$BASE/grelu

mkdir -p $GENOME_O $LABELS $GRELU

echo "Project root: $PROJECT_ROOT"
echo "Species:      $sp"
echo "Genome:       $GENOME"
echo "Annotation:   $ANNOTATION"
echo "Enhancers:    $ENHANCERS"
echo "Labels:       $LABELS"
echo "gReLU output: $GRELU"

##########################################################################################################

# indexing genome assembly

if [[ ! -s $GENOME.fai ]]; then
    samtools faidx $GENOME
fi

# Creating 501 bp centred windows from CRE intervals

if [[ ! -s $LABELS/enhancers_501bp.bed ]]; then

    awk 'BEGIN{OFS="\t"}{
        mid=int(($2+$3)/2);
        print $1,mid,mid+1
    }' \
    $ENHANCERS \
    | bedtools slop -b 250 -g $GENOME.fai -i - | awk 'BEGIN{OFS="\t"} $3-$2 == 501' \
    > $LABELS/enhancers_501bp.bed
fi

# filter N function

filter_N_content () {

    INPUT=$1
    OUTPUT=$2

    awk '
    BEGIN {FS=OFS="\t"}
    NR==1 {
        for (i=1;i<=NF;i++) {
            if ($i ~ /num_N$/) n=i
            if ($i ~ /seq_len$/) len=i
        }
        next
    }
    ($n/$len)<=0.01 {
        print $1,$2,$3
    }
    ' $INPUT > $OUTPUT
}

# calculate gc and N content

if [[ ! -s $LABELS/enhancers.bed ]]; then

    bedtools nuc -fi $GENOME -bed $LABELS/enhancers_501bp.bed \
	    > $LABELS/enhancers_nuc.txt
    filter_N_content $LABELS/enhancers_nuc.txt $LABELS/enhancers.bed

    bedtools nuc -fi $GENOME -bed $LABELS/enhancers.bed \
	    > $LABELS/enhancers_filtered_nuc.txt
fi

# background

if [[ ! -s $LABELS/background.bed ]]; then

    echo "Generating background candidates..."

    bedtools makewindows -g $GENOME.fai -w 501 | awk 'BEGIN{OFS="\t"} $3-$2 == 501' \
    > "$LABELS/genome_windows.bed"

    bedtools intersect -a $LABELS/genome_windows.bed -b $LABELS/enhancers.bed -v \
        > $LABELS/background_candidates.bed

    echo "Filtering background candidates for N content..."

    bedtools nuc -fi $GENOME -bed $LABELS/background_candidates.bed \
    	    > $LABELS/background_candidates_nuc.txt

    filter_N_content $LABELS/background_candidates_nuc.txt $LABELS/background_candidates_filtered.bed

    echo "Recalculating nucleotide composition for filtered candidates..."

    bedtools nuc -fi $GENOME -bed $LABELS/background_candidates_filtered.bed \
        > $LABELS/background_candidates_filtered_nuc.txt

    rm $LABELS/genome_windows.bed
fi

# Creating N-pad matched background set

if [[ ! -s $LABELS/background.bed ]]; then

	$PYTHON $PROJECT_ROOT/src/matchbackground.py \
	    	--enhancers $LABELS/enhancers_filtered_nuc.txt \
	    	--background $LABELS/background_candidates_filtered_nuc.txt \
	    	--output $LABELS/background.bed --max_n_fraction 0.01
fi

# Splitting BED files into training / validation / testing datasets

if [[ ! -s $GRELU/train_enhancers.bed ]]; then

    echo "Splitting $sp BED files into training, validation, and test sets..."

    $PYTHON $PROJECT_ROOT/src/splitbeds.py --species $sp --config $SPLITS --input $LABELS --output $GRELU

else
    echo "Training/validation/test BED files already exist. Skipping split."

fi

##########################################################################################################
# Complete
##########################################################################################################

echo  
echo "============================================================"
echo "$sp label preparation complete"
echo "============================================================"
echo "Labels: $LABELS"
echo "gReLU:  $GRELU"
echo "============================================================"
