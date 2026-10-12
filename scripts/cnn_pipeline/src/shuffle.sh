#!/bin/bash

#SBATCH -J shuffle
#SBATCH --ACCOUNT

#SBATCH --nodes=1
#SBATCH --cpus-per-task=5
#SBATCH --mem=10G
#SBATCH --time=1-00:00:00
#SBATCH --output=./logs/04.slurm-%A_%a.out
#SBATCH --error=./logs/04.slurm-%A_%a.err

set -euo pipefail

CONDA_BIN=$(dirname $(which conda))
CONDA_BASE=$(dirname $CONDA_BIN)

export PATH=$CONDA_BIN:$PATH

source $CONDA_BASE/etc/profile.d/conda.sh

conda activate heliconius_cnn

PROJECT_ROOT=$(pwd)
CONFIG=$PROJECT_ROOT/config/config.yaml

DATABASE=$PROJECT_ROOT/data/databases
GENOME_DIR=$PROJECT_ROOT/$(yq -r '.prediction_species.genome' $CONFIG)

SEQUENCES=$PROJECT_ROOT/data/prediction_species/sequences
ATAC=$PROJECT_ROOT/data/prediction_species/ATAC
CNEE=$PROJECT_ROOT/data/prediction_species/CNEE

RESULTS=$PROJECT_ROOT/results/predictions

module load bedtools2/2.31.1-extr

SP=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $PROJECT_ROOT/config/species_list.txt)
GENOME=$(ls $GENOME_DIR/$SP.assembly.*.fasta)

echo "Running task $SLURM_ARRAY_TASK_ID for $SP"

SHUFFLED=$RESULTS/shuffled/$SP
mkdir -p $SHUFFLED

cut -f1-3 $RESULTS/$SP.predicted_enhancers.merged.bed | grep -v '^chrom' > $SHUFFLED/$SP.predicted_enhancers.merged.bed

for i in $(seq 1 1000); do

        if [[ ! -s $SHUFFLED/Rep.$i.bed ]]; then

                (
                        bedtools shuffle -maxTries 1000000 -chrom -noOverlapping -seed $i -excl $RESULTS/$SP.promoters.1kb.bed \
                                -i $SHUFFLED/$SP.predicted_enhancers.merged.bed -g $GENOME.fai \
                                | bedtools sort -i - > $SHUFFLED/Rep.$i.bed
                ) &
        fi

        if (( i % 5 == 0 )); then
                wait
        fi
done

wait

N=$(ls $SHUFFLED/Rep.*.bed 2>/dev/null | wc -l)

if [[ $N -ne 1000 ]]; then
        echo "ERROR: $SP has only $N/1000 shuffled replicates"
        exit 1
fi

ls -1 $SHUFFLED/Rep.*.bed | sort -V > $SHUFFLED/$SP.enhancers.shuffled.fofn

echo "$SP complete"
