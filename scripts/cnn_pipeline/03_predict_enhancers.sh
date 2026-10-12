#!/bin/bash

#SBATCH -J grelu
#SBATCH --ACCOUNT

#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --mem=10G
#SBATCH --time=3-00:00:00
#SBATCH --output=./logs/03.slurm-%j.out
#SBATCH --error=./logs/03.slurm-%j.err

##########################################################################################################
############################ Predict enhancer regions using gReLU CNN model ##############################
##########################################################################################################

# Inputs:
#   - Trained gReLU CNN model
#   - Prediction species genome FASTA assemblies
#   - Prediction species ATAC-seq peak BED files
#
# Outputs:
#   - Genome-wide fixed-width prediction windows
#   - Predicted enhancer BED files and prediction scores
#   - Predicted enhancer overlaps with ATAC-seq peaks
#   - Permutation test results assessing enrichment of predictions within ATAC peaks

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

##########################################################################################################
# gReLU setup
##########################################################################################################

echo
echo "============================================================"
echo "gReLU pipeline"
echo "============================================================"

PROJECT_ROOT=$(pwd)
CONFIG=$PROJECT_ROOT/config/config.yaml

GENOME_DIR=$PROJECT_ROOT/$(yq -r '.prediction_species.genome' $CONFIG)

SEQUENCES=$PROJECT_ROOT/data/prediction_species/sequences
ATAC=$PROJECT_ROOT/data/prediction_species/ATAC
ANNOTATION=$PROJECT_ROOT/data/prediction_species/annotation

RESULTS=$PROJECT_ROOT/results/predictions
REGMAP=$RESULTS/regmap

mkdir -p $SEQUENCES $ATAC $RESULTS $REGMAP

echo "Project root: $PROJECT_ROOT"
echo "Genome dir:   $GENOME_DIR"
echo "Sequences:    $SEQUENCES"
echo "ATAC:         $ATAC"
echo "Results:      $RESULTS"
echo "RegMap:       $REGMAP"

##########################################################################################################

echo "Preparing query genomes..."

for GENOME in $GENOME_DIR/*.fasta; do

	SP=$(basename $GENOME | cut -d'.' -f1)
	echo "Processing $SP..."

	if [[ ! -s $GENOME.fai ]]; then
		samtools faidx $GENOME
	fi

	if [[ ! -s $SEQUENCES/$SP.501bp_windows.bed ]]; then

		bedtools makewindows -g $GENOME.fai -w 501 -s 100 \
			> $SEQUENCES/$SP.501bp_windows.bed
	fi
	echo "Finished creating genomic windows for $SP."
done

echo "Query genomes prepared."

MODEL=$PROJECT_ROOT/models/grelu/best.ckpt

for GENOME in $GENOME_DIR/*.fasta; do
        SP=$(basename $GENOME | cut -d'.' -f1)

	if [[ ! -s $RESULTS/$SP.predicted_enhancers.bed ]]; then

		mkdir -p $RESULTS
		echo "Running enhancer prediction model..."

		$PYTHON $PROJECT_ROOT/src/predict_grelu.py --fasta $GENOME --bed $SEQUENCES/$SP.501bp_windows.bed \
			--model $MODEL --output $RESULTS/$SP.predicted_enhancers.bed
	fi

	if [[ ! -s $RESULTS/$SP.predicted_enhancers.merged.bed ]]; then

                echo "Removing ±1 kb TSS regions and merging predicted enhancer windows for $SP..."

                bedtools slop \
                        -i $ANNOTATION/$SP.annotation.TSS.bed -g $GENOME.fai -b 1000 > $RESULTS/$SP.promoters.1kb.bed
                bedtools intersect -v -a $RESULTS/$SP.predicted_enhancers.bed -b $RESULTS/$SP.promoters.1kb.bed | \
                        sort -k1,1 -k2,2n | \
                        bedtools merge -i - -c 4,4 -o max,count > $RESULTS/$SP.predicted_enhancers.merged.bed
        fi
done


echo -e "species\toriginal_merged\tfiltered_merged\toriginal_bp\tfiltered_bp\tbp_removed\tpercent_bp_removed" \
        > $RESULTS/promoter_filter_merged_summary.tsv

for GENOME in $GENOME_DIR/*.fasta; do

        SP=$(basename $GENOME | cut -d'.' -f1)

        ORIGINAL=$(sort -k1,1 -k2,2n $RESULTS/$SP.predicted_enhancers.bed | \
                bedtools merge -i - | wc -l)

        FILTERED=$(wc -l < $RESULTS/$SP.predicted_enhancers.merged.bed)

        ORIGINAL_BP=$(sort -k1,1 -k2,2n $RESULTS/$SP.predicted_enhancers.bed | \
                bedtools merge -i - | \
                awk '{sum += $3-$2} END {print sum}')

        FILTERED_BP=$(awk '{sum += $3-$2} END {print sum}' \
                $RESULTS/$SP.predicted_enhancers.merged.bed)

        REMOVED_BP=$((ORIGINAL_BP - FILTERED_BP))

        PERCENT_BP=$(awk -v r=$REMOVED_BP -v o=$ORIGINAL_BP \
                'BEGIN {printf "%.2f", (r/o)*100}')

        echo -e "$SP\t$ORIGINAL\t$FILTERED\t$ORIGINAL_BP\t$FILTERED_BP\t$REMOVED_BP\t$PERCENT_BP" \
                >> $RESULTS/promoter_filter_merged_summary.tsv
done

########### Determining overlap of predicted enhancers with ATAC peaks ############

for SP in Hcha Hdem Hmel; do

	ATAC_overlap=$RESULTS/ATAC_overlap
	mkdir -p $ATAC_overlap

  	SUMMARY=$ATAC_overlap/$SP.ATAC_enhancer_overlap_summary.tsv
        echo -e "species\tcondition\ttotal_ATAC_peaks\toverlap_peaks\tpercent_overlap" > $SUMMARY

	PEAKS=($ATAC/$SP.*.bestMACS3.narrowPeak)

	if [[ ! -s $ATAC/$SP.ATAC.merged.bed ]]; then
		
		echo "Creating merged ATAC set for $SP..."

		cat "${PEAKS[@]}" | cut -f1-3 | sort -k1,1 -k2,2n | bedtools merge -i - \
			> $ATAC/$SP.ATAC.merged.bed
	fi

	for PEAK in ${PEAKS[@]}; do

		BASE=$(basename "$PEAK")
                Cond=${BASE#${SP}.}
                Cond=${Cond%.bestMACS3.narrowPeak}

		if [[ ! -s $ATAC_overlap/$SP.enhancer_${Cond}_overlap.bed ]]; then

			echo "Overlapping condition $Cond with predicted $SP enhancers..."

			bedtools intersect -u -f 0.5 -a $PEAK -b $RESULTS/$SP.predicted_enhancers.merged.bed \
				> $ATAC_overlap/$SP.enhancer_${Cond}_overlap.bed

			TOTAL=$(wc -l < $PEAK)
			OVERLAP=$(wc -l < $ATAC_overlap/$SP.enhancer_${Cond}_overlap.bed)
			PERCENT=$(awk -v o=$OVERLAP -v t=$TOTAL 'BEGIN {printf "%.2f", (o/t)*100}')

			echo -e "$SP\t$Cond\t$TOTAL\t$OVERLAP\t$PERCENT" >> $SUMMARY
		fi
	done

	if [[ ! -s $ATAC_overlap/$SP.enhancer_all_overlap.bed ]]; then

		echo "Overlapping merged ATAC conditions with predicted $SP enhancers..."

		bedtools intersect -u -f 0.5 -a $ATAC/$SP.ATAC.merged.bed -b $RESULTS/$SP.predicted_enhancers.merged.bed \
			> $ATAC_overlap/$SP.enhancer_all_overlap.bed

		TOTAL=$(wc -l < $ATAC/$SP.ATAC.merged.bed)
		OVERLAP=$(wc -l < $ATAC_overlap/$SP.enhancer_all_overlap.bed)
		PERCENT=$(awk -v o=$OVERLAP -v t=$TOTAL 'BEGIN {printf "%.2f", (o/t)*100}')

		echo -e "$SP\tAllMerged\t$TOTAL\t$OVERLAP\t$PERCENT" >> $SUMMARY
	fi

	### Permutation test ####
        
        # create shuffled enhancer background

        SHUFFLED=$RESULTS/shuffled/$SP
        mkdir -p $SHUFFLED

        if [[ ! -s $SHUFFLED/$SP.enhancers.shuffled.fofn ]]; then

                cut -f1-3 $RESULTS/$SP.predicted_enhancers.merged.bed | grep -v '^chrom' > $SHUFFLED/$SP.predicted_enhancers.merged.bed

                for i in $(seq 1 1000); do

                        bedtools shuffle -maxTries 1000000 -chrom -noOverlapping -seed $i -excl $RESULTS/$SP.promoters.1kb.bed \
                                -i $SHUFFLED/$SP.predicted_enhancers.merged.bed -g $GENOME_DIR/$SP.assembly.*.fasta.fai \
                                | bedtools sort -i - > $SHUFFLED/Rep.$i.bed
                done

                ls -1 $SHUFFLED/Rep.*.bed | sort -V > $SHUFFLED/$SP.enhancers.shuffled.fofn
        fi

        RUNPERMUTATION=1

        if [[ $RUNPERMUTATION == 1 ]]; then

                echo "Running permutation test for $SP..."

                OBS=$(bedtools intersect -u -f 0.5 -a $ATAC/$SP.ATAC.merged.bed -b $RESULTS/$SP.predicted_enhancers.merged.bed | wc -l)

                echo -e "rep\toverlap" > $ATAC_overlap/$SP.ATAC_permutations.tsv

                rep=0
                while read SHUFFLE; do
                        ((++rep))

                        COUNT=$(bedtools intersect -u -f 0.5 -a $ATAC/$SP.ATAC.merged.bed -b $SHUFFLE | wc -l)

                        echo -e "$rep\t$COUNT" >> $ATAC_overlap/$SP.ATAC_permutations.tsv

                        if (( rep % 100 == 0 )); then
                                echo "Completed $rep permutations..."
                        fi

                done < $SHUFFLED/$SP.enhancers.shuffled.fofn

                MEAN=$(awk 'NR>1 {sum+=$2; n++} END {print sum/n}' $ATAC_overlap/$SP.ATAC_permutations.tsv)
                FOLD=$(awk -v obs=$OBS -v mean=$MEAN 'BEGIN {print obs/mean}')
                P=$(awk -v obs=$OBS '
                        NR>1 {
                                n++
                                if ($2 >= obs) x++
                        }
                        END {
                                print (x+1)/(n+1)
                        }' $ATAC_overlap/$SP.ATAC_permutations.tsv)

                echo -e "species\tobserved\tpermutation_mean\tfold_enrichment\tempirical_p" \
                        > $ATAC_overlap/$SP.ATAC_permutation_summary.tsv

                echo -e "$SP\t$OBS\t$MEAN\t$FOLD\t$P" \
                        >> $ATAC_overlap/$SP.ATAC_permutation_summary.tsv

                echo "Observed: $OBS"
                echo "Permutation mean: $MEAN"
                echo "Fold enrichment: $FOLD"
                echo "Empirical p: $P"
	fi
done
