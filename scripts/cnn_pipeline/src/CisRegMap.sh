#!/bin/bash

#SBATCH -J CisRegMap.sh
#SBATCH --ACCOUNT
#SBATCH --nodes=1
#SBATCH --cpus-per-task=5
#SBATCH --mem=64G
#SBATCH --time=5-06:00:00
#SBATCH --output=./logs/regmap_logs/tomtom.slurm-%j.out
#SBATCH --error=./logs/regmap_logs/tomtom.slurm-%j.err

##################################################################################
########## CisRegMap and motif analysis for CNN enhancer predictions #############
##################################################################################

set -eo pipefail

source /user/work/wd24922/miniforge3/etc/profile.d/conda.sh
conda activate CisRegMapEnv

set -u

export PATH="$CONDA_PREFIX/bin:$PATH"
hash -r

PYTHON="$CONDA_PREFIX/bin/python"

echo "Conda environment: $CONDA_PREFIX"
echo "Python: $(which python)"
python --version

module load languages/R/4.5.1

########## CisRegMap.sh ##########

PROJECT_ROOT=$(pwd)

CONFIG=$PROJECT_ROOT/config/config.yaml
GENOME_DIR=$PROJECT_ROOT/data/prediction_species/genomes
ANNOTATION=$PROJECT_ROOT/data/prediction_species/annotation
RESULTS=$PROJECT_ROOT/results/predictions
RSA=$PROJECT_ROOT/src/matrix-clustering_stand-alone
MEME=$PROJECT_ROOT/src/meme-5.5.6
DB=$PROJECT_ROOT/data/databases/JASPAR2024_CORE_insects_non-redundant_pfms_meme.txt

export PATH=$MEME/bin:$PATH

echo "MEME: $MEME"
which tomtom
tomtom -version

SP=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$PROJECT_ROOT/config/species_list.txt")

echo "Processing $SP"

$PROJECT_ROOT/src/get_TSSandGeneBody.py -i $ANNOTATION/$SP.annotation.CAT.gff3 -o $ANNOTATION/$SP.annotation.TSS.bed

REGMAP=$RESULTS/regmap
mkdir -p $REGMAP

SHUFFLED=$REGMAP/shuffled/$SP
mkdir -p $SHUFFLED

if [[ ! -s $REGMAP/shuffled/$SP.enhancers.shuffled.fofn ]]; then

	cut -f1-3 $RESULTS/$SP.predicted_enhancers.merged.bed \
    		| grep -v '^chrom' \
    		> $SHUFFLED/$SP.predicted_enhancers.merged.bed

	for i in $(seq 1 1000); do
    		bedtools shuffle \
			-maxTries 1000000 \
			-f 0.95 \
			-chrom \
			-i $SHUFFLED/$SP.predicted_enhancers.merged.bed \
			-g $GENOME_DIR/$SP.assembly.*.fasta.fai \
	    		| bedtools sort -i - \
	    		> $SHUFFLED/Rep.$i.bed
	done

	ls -1 $SHUFFLED/Rep.*.bed | sort -V \
    		> $REGMAP/shuffled/$SP.enhancers.shuffled.fofn
fi

OUTPUT=$REGMAP/$SP

if [[ ! -s $OUTPUT/Class.tsv ]]; then

	mkdir -p $OUTPUT/work
	cd $OUTPUT/work

	$PYTHON -u $PROJECT_ROOT/src/CisRegMap.py \
	    	-t $ANNOTATION/$SP.annotation.TSS.bed \
		-e $ANNOTATION/$SP.annotation.exons.bed \
		-i $ANNOTATION/$SP.annotation.introns.bed \
	    	-p $RESULTS/$SP.predicted_enhancers.merged.bed \
	    	-r $REGMAP/shuffled/$SP.enhancers.shuffled.fofn
	mv stacked_peak_distribution.pdf $OUTPUT/stacked_peak_distribution.pdf
	mv AssGene.tsv $OUTPUT/AssGene.tsv
	mv Class.tsv $OUTPUT/Class.tsv
fi

echo "$SP regulatory mapping complete"


########## Motif analysis ##########

MOTIF=$RESULTS/motifs
mkdir -p $MOTIF

echo "Conducting motif analysis for $SP..."

if [[ ! -s $MOTIF/$SP.HomerEnrichedMotif/homerResults.html ]]; then

	findMotifsGenome.pl $RESULTS/$SP.predicted_enhancers.merged.bed \
		$GENOME_DIR/$SP.assembly.*.fasta \
		$MOTIF/$SP.HomerEnrichedMotif \
		-bg $SHUFFLED/Rep.1.bed \
		-mset insects \
		-size given
fi

if [[ ! -s $MOTIF/$SP.HomerEnrichedMotif/unified.motifs ]]; then

	cat $MOTIF/$SP.HomerEnrichedMotif/homerMotifs.all.motifs $MOTIF/$SP.HomerEnrichedMotif/knownResults/known*.motif >> $MOTIF/$SP.HomerEnrichedMotif/unified.motifs
fi

if [[ ! -s $MOTIF/$SP.cluster_tables/summary_table.tab ]]; then

	Rscript $RSA/matrix-clustering.R \
		--M $MOTIF/$SP.HomerEnrichedMotif/unified.motifs \
		--C $SP \
		--F homer \
		-o $MOTIF/$SP.cluster \
		--export_heatmap TRUE
fi

mkdir -p $MOTIF/$SP.tomtom

if [[ ! -s $MOTIF/$SP.tomtom/${SP}_output ]]; then

	if [[ ! -s $MOTIF/$SP.tomtom/$SP.unified.meme ]]; then
		
		$PYTHON $PROJECT_ROOT/src/homer2meme.py \
			$MOTIF/$SP.HomerEnrichedMotif/unified.motifs \
			$MOTIF/$SP.tomtom/$SP.unified.meme
	fi

	echo "Running Tomtom for $SP..."

	$MEME/bin/tomtom \
		-no-ssc \
		-oc $MOTIF/$SP.tomtom/${SP}_output \
		-min-overlap 5 \
		-dist pearson \
		-evalue \
		-thresh 0.1 \
		$MOTIF/$SP.tomtom/$SP.unified.meme \
		$DB
fi
 
mkdir -p $MOTIF/families

if [[ ! -s $MOTIF/families/${SP}_families.tsv ]]; then

	echo "Concatenating tomtom.tsv with annotation metadata for $SP..."
	
	$PYTHON $PROJECT_ROOT/src/annotate_tomtom.py \
	    	-m $PROJECT_ROOT/data/databases/JASPAR2024_TF_metadata.tsv \
	    	-i $MOTIF/$SP.tomtom/${SP}_output/ \
	    	-o $MOTIF/families/${SP}_families.tsv

fi
