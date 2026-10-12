#!/bin/bash

#SBATCH -J CNEE
#SBATCH --ACCOUNT

#SBATCH --nodes=1
#SBATCH --cpus-per-task=5
#SBATCH --mem=10G
#SBATCH --time=14-00:00:00
#SBATCH --output=./logs/04.slurm-%j.out
#SBATCH --error=./logs/04.slurm-%j.err

##########################################################################################################
########################## Intersect predicted enhancers with conserved CNEEs ############################
##########################################################################################################

# Inputs:
#   - Predicted enhancer BED files for prediction species
#   - H. melpomene CNEE reference BED file
#   - Multi-species HAL genome alignment
#   - Prediction species genome assemblies
#
# Outputs:
#   - H. melpomene CNEEs lifted to each prediction species
#   - Predicted enhancer overlaps with lifted CNEEs
#   - Species-specific CNEE-associated predicted enhancer BED files

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

HAL=/user/work/wd24922/software/hal
export PATH=$HAL/bin:$PATH

echo
echo "============================================================"
echo "Environment"
echo "============================================================"

echo "CONDA_EXE=$CONDA_EXE"
echo "CONDA_PREFIX=$CONDA_PREFIX"
echo "Python=$PYTHON"
echo "Node=$(hostname)"
echo "SLURM_JOB_ID=${SLURM_JOB_ID:-unset}"
echo "CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES:-unset}"

echo
echo "============================================================"
echo "HAL tools"
echo "============================================================"

echo "halLiftover: $(which halLiftover)"

##########################################################################################################
# Project paths
##########################################################################################################

echo
echo "============================================================"
echo "CNEE conservation pipeline"
echo "============================================================"

PROJECT_ROOT=$(pwd)
CONFIG=$PROJECT_ROOT/config/config.yaml

DATABASE=$PROJECT_ROOT/data/databases
GENOME_DIR=$PROJECT_ROOT/$(yq -r '.prediction_species.genome' $CONFIG)

SEQUENCES=$PROJECT_ROOT/data/prediction_species/sequences
ATAC=$PROJECT_ROOT/data/prediction_species/ATAC
CNEE=$PROJECT_ROOT/data/prediction_species/CNEE

RESULTS=$PROJECT_ROOT/results/predictions

mkdir -p $CNEE

echo "Project root: $PROJECT_ROOT"
echo "Genome dir:   $GENOME_DIR"
echo "Database:     $DATABASE"
echo "CNEE:         $CNEE"
echo "Predictions:  $RESULTS"

##########################################################################################################

# CNEE HalLitover

CACTUS=$DATABASE/63Nymphalidae.Cactus-v1.5.hal

CNEE_overlap=$RESULTS/CNEE_overlap
mkdir -p $CNEE_overlap

if [[ ! -s $CNEE/Hmel.CNEEs.bed ]]; then

        echo "Preparing Hmel CNEEs..."

        cut -f1-4 $CNEE/Hmel.phyloacc.allCNEEs.bed > $CNEE/Hmel.CNEEs.bed
fi

echo -e "species\ttotal_CNEEs\tCNEEs_overlapping_enhancers\tpercent_CNEEs_overlap\ttotal_enhancers\tenhancers_overlapping_CNEEs\tpercent_enhancers_overlap" \
        > $CNEE_overlap/CNEE_enhancer_overlap_summary.tsv

for GENOME in $GENOME_DIR/*.fasta; do

        SP=$(basename $GENOME | cut -d'.' -f1)
        if [[ ! -s $CNEE/$SP.CNEEs.bed ]]; then

                echo "Lifting Hmel CNEEs to $SP..."
                halLiftover --noDupes $CACTUS Hmel $CNEE/Hmel.CNEEs.bed $SP $CNEE/$SP.CNEEs.bed
        fi

        if [[ ! -s $CNEE_overlap/$SP.CNEEs_overlapping_enhancers.bed ]]; then

                echo "Intersecting CNEEs with predicted $SP enhancers..."
                bedtools intersect -u -a $CNEE/$SP.CNEEs.bed -b $RESULTS/$SP.predicted_enhancers.merged.bed \
                        > $CNEE_overlap/$SP.CNEEs_overlapping_enhancers.bed
        fi

        if [[ ! -s $CNEE_overlap/$SP.enhancers_overlapping_CNEEs.bed ]]; then

                echo "Intersecting predicted $SP enhancers with CNEEs..."
                bedtools intersect -u -a $RESULTS/$SP.predicted_enhancers.merged.bed -b $CNEE/$SP.CNEEs.bed \
                        > $CNEE_overlap/$SP.enhancers_overlapping_CNEEs.bed
        fi

        TOTAL_CNEE=$(cut -f4 $CNEE/$SP.CNEEs.bed | sort -u | wc -l)
        CNEE_OVERLAP=$(cut -f4 $CNEE_overlap/$SP.CNEEs_overlapping_enhancers.bed | sort -u | wc -l)
        CNEE_PERCENT=$(awk -v o=$CNEE_OVERLAP -v t=$TOTAL_CNEE 'BEGIN {printf "%.2f", (o/t)*100}')

        TOTAL_ENH=$(wc -l < $RESULTS/$SP.predicted_enhancers.merged.bed)
        ENH_OVERLAP=$(wc -l < $CNEE_overlap/$SP.enhancers_overlapping_CNEEs.bed)
        ENH_PERCENT=$(awk -v o=$ENH_OVERLAP -v t=$TOTAL_ENH 'BEGIN {printf "%.2f", (o/t)*100}')

        echo -e "$SP\t$TOTAL_CNEE\t$CNEE_OVERLAP\t$CNEE_PERCENT\t$TOTAL_ENH\t$ENH_OVERLAP\t$ENH_PERCENT" \
                >> $CNEE_overlap/CNEE_enhancer_overlap_summary.tsv
done
	#### permuations ####

echo -e "species\tdirection\tobserved\tpermutation_mean\tfold_enrichment\tempirical_p" > $CNEE_overlap/CNEE_permutation_summary.tsv

for GENOME in $GENOME_DIR/*.fasta; do

	SP=$(basename $GENOME | cut -d'.' -f1)
	SHUFFLED=$RESULTS/shuffled/$SP
	mkdir -p $SHUFFLED

	if [[ ! -s $SHUFFLED/$SP.enhancers.shuffled.fofn ]]; then

		cut -f1-3 $RESULTS/$SP.predicted_enhancers.merged.bed | grep -v '^chrom' > $SHUFFLED/$SP.predicted_enhancers.merged.bed

		for i in $(seq 1 1000); do

			bedtools shuffle -maxTries 1000000 -chrom -seed $i -excl $RESULTS/$SP.promoters.1kb.bed \
				-i $SHUFFLED/$SP.predicted_enhancers.merged.bed -g $GENOME_DIR/$SP.assembly.*.fasta.fai \
				| bedtools sort -i - > $SHUFFLED/Rep.$i.bed
		done
		
		ls -1 $SHUFFLED/Rep.*.bed | sort -V > $SHUFFLED/$SP.enhancers.shuffled.fofn
	fi

	RUNPERMUTATION=1

 	if [[ $RUNPERMUTATION == 1 ]]; then

                echo "Running permutation test for $SP..."

                OBS_CNEE=$(bedtools intersect -wa -u -a $CNEE/$SP.CNEEs.bed -b $RESULTS/$SP.predicted_enhancers.merged.bed | cut -f4 | sort -u | wc -l)
                OBS_ENH=$(bedtools intersect -u -a $RESULTS/$SP.predicted_enhancers.merged.bed -b $CNEE/$SP.CNEEs.bed | wc -l)

                if [[ ! -s $CNEE_overlap/$SP.CNEE_permutations.tsv || $(wc -l < $CNEE_overlap/$SP.CNEE_permutations.tsv) -ne 1001 ]]; then

                        echo -e "rep\tCNEE_overlap\tenhancer_overlap" > $CNEE_overlap/$SP.CNEE_permutations.tsv

                        rep=0

                        while read SHUFFLE; do

                                ((++rep))

                                COUNT_CNEE=$(bedtools intersect -wa -u -a $CNEE/$SP.CNEEs.bed -b $SHUFFLE | cut -f4 | sort -u | wc -l)
                                COUNT_ENH=$(bedtools intersect -u -a $SHUFFLE -b $CNEE/$SP.CNEEs.bed | wc -l)

                                echo -e "$rep\t$COUNT_CNEE\t$COUNT_ENH" >> $CNEE_overlap/$SP.CNEE_permutations.tsv

                                if (( rep % 100 == 0 )); then
                                        echo "Completed $rep permutations..."
                                fi

                        done < $SHUFFLED/$SP.enhancers.shuffled.fofn
                else
                        echo "$SP permutation file already complete. Skipping permutations."
                fi

                MEAN_CNEE=$(awk 'NR>1 {sum+=$2; n++} END {print sum/n}' $CNEE_overlap/$SP.CNEE_permutations.tsv)
                FOLD_CNEE=$(awk -v obs=$OBS_CNEE -v mean=$MEAN_CNEE 'BEGIN {print obs/mean}')
                P_CNEE=$(awk -v obs=$OBS_CNEE '
                        NR>1 {
                                n++
                                if ($2 >= obs) x++
                        }
                        END {
                                print (x+1)/(n+1)
                        }' $CNEE_overlap/$SP.CNEE_permutations.tsv)

                MEAN_ENH=$(awk 'NR>1 {sum+=$3; n++} END {print sum/n}' $CNEE_overlap/$SP.CNEE_permutations.tsv)
                FOLD_ENH=$(awk -v obs=$OBS_ENH -v mean=$MEAN_ENH 'BEGIN {print obs/mean}')
                P_ENH=$(awk -v obs=$OBS_ENH '
                        NR>1 {
                                n++
                                if ($3 >= obs) x++
                        }
                        END {
                                print (x+1)/(n+1)
                        }' $CNEE_overlap/$SP.CNEE_permutations.tsv)

                echo -e "$SP\tCNEE_to_enhancer\t$OBS_CNEE\t$MEAN_CNEE\t$FOLD_CNEE\t$P_CNEE" >> $CNEE_overlap/CNEE_permutation_summary.tsv
                echo -e "$SP\tenhancer_to_CNEE\t$OBS_ENH\t$MEAN_ENH\t$FOLD_ENH\t$P_ENH" >> $CNEE_overlap/CNEE_permutation_summary.tsv
        fi
done
