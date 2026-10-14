#!/bin/bash

#### key paths ####

path="/user/work/wd24922/Project/Practice"
in_path="$path/Input"
out_path="$path/Output"
assembly_path="$path/Assemblies"
hal_path="/user/home/wd24922/work/Project/AssembledGenomes"

#### software ####

APPS="/user/work/wd24922/miniforge3/bin"
APPS2="/user/work/tk19812/software"
SOFTWARE="/user/work/wd24922/software"
PICARD_HOME="$APPS2/PicardTools"
BIOSCRIPT="$APPS2/bioscript"
BIOSCRIPT2="$SOFTWARE/bioscript2"
BIOSCRIPT3="/user/home/tk19812/software/bioscript"
DEEPTOOLS="$APPS2/deepTools/deeptools"
MINICONDA_PYTHON="/user/home/wd24922/miniconda3/bin/python"
HOMER="$APPS2/Homer/bin"
RSA="$SOFTWARE/matrix-clustering_stand-alone"
MEME2="$APPS2/meme-5.5.6/bin"
CHROM="$SOFTWARE/chrombpnet/chrombpnet"
MOTIF="$path/databases"

THREADS=20
TRIMMER=trim_galore
MAPPER=bwa

module load bwa
module load sratoolkit
module load samtools
module load bedtools2
module load languages/perl/5.38.2
module load languages/R/4.4.3

# Initialise conda for non-interactive shells
source /user/work/wd24922/miniforge3/etc/profile.d/conda.sh
conda activate CisRegMapEnv

export PATH="/user/work/wd24922/miniforge3/envs/CisRegMapEnv/bin:$PATH"
export PATH="$HOMER:$PATH"
export PATH="$APPS:$PATH"
MEME=/user/work/wd24922/software/meme-5.5.9
export PATH="$MEME/bin:$PATH"

set -eo pipefail

################################################# Characterisation of cis-regulatory elements across spatiotemporal conditions pipeline ########################################################

# Controls (0 = off / 1 = on)

#### data preprocessing ####

DOWNLOADFASTQ=0 # Downloads data from SRR ids
RUNBUSCO=0 # Scores genome confidence using compleasm and busco
RUNPREPROCESSING=0  # Performs data preprocessing and mapping
RUNCHECK=0 # Quality checks replicate similarity
RUNFRAGMENT=0 # Performs ATAC seq fragment analysis across individual conditions and universal datasets
RUNPEAKCALLING=0 # Performs peak calling with MACS3 on all condition types
RUNQUALITYCHECK=0  # Quality checking of peaks with deeptools

#### characterising cis regulatory elements ####

RUNCONSENSUS=0 # Creates consensus ATAC, ac, and met dataset for creation of promoter and enhancer files
RUNINTERSECT=0 # Creates promoter and enhancer bed files using bedtools intersect for downstream analysis alongside permutation tests for statistical significance of overlap
RUNCISCHECK=0 # Permutations tests and plot profile
RUNREGMAP=0 # Runs cis regulatory mapping of cis regulatory elements and permutations test (for slurm array when running cisregmap, submit with N=$(wc -l < /user/work/wd24922/Project/Practice/Output/CisData/regmap/regmap_tasks.tsv) sbatch --array=1-$N%6 cis_characterisation.sh)
RUNGENEANNOTATION=0 # Determines if putative enhancers and promoters are associated with known Heliconius wing development genes using orthologs from melpomene

#### characterising transcription factors ####

RUNMOTIF=0 # Performs motif identification of de novo and known motifs using HOMER on promoter and enhancer bed files
RUNGC=0 # Running bedtools nuc to assess GC content of cis-regulatory beds compared to genome
RUNCLUSTER=0 # Clustering of homer motifs via R-matrix clustering to assess possible TF families
RUNTOBIAS=0 # Takes unified motif from homer and runs with TOBIAS to correct Tn5 insertion bias, perform footprinting to identify active motifs, and further performs bindetect to identify differential expression of promoter and enhancer transcription factors
RUNTFMATCH=0 # Passes motifs through TOMTOM to identify closest associated transcription factors using JASPAR insect database

#### deseq2 differential expression analysis ####

RUNDESEQ2=0 # Creates featureCounts matrix, column sample data, and runs deseq2 for comparison of ATAC, ac, and met promoter and enhancers marks across spatiotemporal conditions

###########################################################################################################################################################################################

# submit job for each species

mapfile -t SPECIES < <(grep -vw "Species" "$in_path/D3SeqData.tsv" | cut -f1 | sort -u)

if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    echo "SLURM_ARRAY_TASK_ID not set. Running locally for species index 0."
    task_id=0
else
    task_id=$SLURM_ARRAY_TASK_ID
fi

sp=${SPECIES[$task_id]}

echo "Running pipeline for $sp"

# ./run_pipeline.sh "$sp" to test code with individual species


#######################################################################################################################################################################################
####################################################################### ATAC-seq and ChiP-seq data preprocessing ######################################################################
#######################################################################################################################################################################################

############################################ Downloading Fastq data from SRR ids ####################################################

# using fasterq-dump for obtaining SRR ids

if [[ $DOWNLOADFASTQ == 1 ]]; then

    	for dataset in D3 OM; do

		fastqData="$in_path/FastqData/${dataset}FastqData"
		mkdir -p "$fastqData"

		mapfile -t SRR_LIST < <(
    		awk -F'\t' -v sp=$sp 'NR>1 && $1==sp {print $5}' $in_path/${dataset}SeqData.tsv | sort -u
		)

		for SRR in ${SRR_LIST[@]}; do

	    		fq1=$fastqData/${SRR}_1.fastq
	    		fq2=$fastqData/${SRR}_2.fastq

	    		if [[ -s $fq1 && -s $fq2 ]]; then
				echo "FASTQs for $SRR already exist, skipping"
				continue
	    		fi

	    		echo "Downloading $SRR for species $sp..."

	    		prefetch $SRR

	    		fasterq-dump -e $THREADS -p --split-files $SRR -O $fastqData
		done
    	done
fi

########################################### Scoring assembled genome confidence ################################################

if [[ $RUNBUSCO == 1 ]]; then

	ASSEMBLY=$(ls -1 $assembly_path/$sp.assembly.*.fasta | head -n 1)

	if [[ ! -f $ASSEMBLY ]]; then
		echo "Assembly not found for $sp, skipping"
		continue
	fi

	number=$(cat $ASSEMBLY | grep -v ">" | tr -d '\n' | wc -c)
	echo $number

	export BUSCO_DOWNLOADS=$out_path/genomes/busco

	BUSCO_OUT=$out_path/genomes/busco/$sp
	mkdir -p $BUSCO_OUT

	if [[ ! -s $BUSCO_OUT ]]; then

		busco -i $ASSEMBLY -l lepidoptera_odb12 -o $BUSCO_OUT -m genome -c $THREADS -f --tar --download_path $BUSCO_DOWNLOADS
	fi

	mkdir -p $out_path/genomes/contig

	if [[ ! -s $out_path/genomes/contig/$sp.txt ]]; then

		$APPS/assemblyStatistics $ASSEMBLY 
	fi
fi

################################## Data preprocessing mapping ChIP-seq and ATAC-Seq Data against indexed genomes #############################

# using fastqc, trimgalore, bwa, picardtools, and samtools 

if [[ $RUNPREPROCESSING == 1 ]]; then   	

	for dataset in D3 OM; do	
		
		Libs=$(grep -w $sp $in_path/${dataset}SeqData.tsv | cut -f2 | sort -V ) # finds unique library identifiers associated with each species

		if [[ -z $Libs ]]; then
		    echo "No libraries found for species $sp"
		    continue
		fi
	
		for lib in $Libs; do # first sections iterates over each library and extracts the tissue, condition, and srrid for each library
			Tissue=$(grep -w $lib $in_path/${dataset}SeqData.tsv | cut -f3 )
			Cond=$(grep -w $lib $in_path/${dataset}SeqData.tsv | cut -f4 )
			SRR=$(awk -F'\t' -v lib="$lib" '$2==lib {print $5}' "$in_path/${dataset}SeqData.tsv" | tr -d '\r')
			rep=$(awk -F'\t' -v lib="$lib" '$2==lib {print $6}' "$in_path/${dataset}SeqData.tsv" | head -n1)	

			# defining file paths for inputs
			ReadA=$in_path/FastqData/${dataset}FastqData/${SRR}_1.fastq
			ReadB=$in_path/FastqData/${dataset}FastqData/${SRR}_2.fastq
			
			if [[ ! -f $ReadA ]]; then
				echo "Input file $ReadA not found for species $sp, library $lib"
				continue
			fi

			if [[ ! -f $ReadB ]]; then
				echo "Input file $ReadB not found for species $sp, library $lib"
				continue
			fi
	
			cons=$out_path/conservation
			preprocess=$cons/preprocessing
			OUTPREX=$dataset.$sp.$lib.$Tissue.$Cond
	
			mkdir -p $cons
			mkdir -p $preprocess
			
			source $(conda info --base)/etc/profile.d/conda.sh
			conda activate cutadapt
			module add fastqc/0.12.1 # adds fastqc

			if [[ ! -s $preprocess/${OUTPREX}_val_1.fq.gz ]]; then

				echo "Running trim galore on $OUTPREX..."
				
				$APPS2/TrimGalore-0.6.10/trim_galore --paired --cores 8 --gzip --length 20 --trim-n --basename $(basename $OUTPREX) -o $preprocess $ReadA $ReadB #runs trim galore
			fi

			trimmedA=$preprocess/${OUTPREX}_val_1.fq.gz	
			trimmedB=$preprocess/${OUTPREX}_val_2.fq.gz

			if [[ ! -f $trimmedA || ! -f $trimmedB ]]; then
				echo "Trimmed files not found"
				continue
			fi
			
				echo "Running FastQC on trimmed reads..."

				mkdir -p $preprocess/fastqc
		
				fastqc -t $THREADS -o $preprocess/fastqc $trimmedA $trimmedB # runs fastqc
						
			multiqc $preprocess/fastqc/.

			if [[ ! -s $preprocess/$OUTPREX.$TRIMMER.$MAPPER.Aligned.bam ]]; then

				echo "Running bwa alignment on $OUTPREX..."

				INDEX=`ls -1 $assembly_path/$sp.assembly.*.fasta`

				bwa mem -t $THREADS $INDEX $trimmedA $trimmedB \
					| samtools view -@ $THREADS -b -S \
					| samtools sort -o $preprocess/$OUTPREX.$TRIMMER.$MAPPER.Aligned.bam # aligns reads with bwa and converts to bam format
			fi

			BAM=$preprocess/$OUTPREX.$TRIMMER.$MAPPER.Aligned.bam

			if [[ ! -s ${BAM%.bam}.filtered.bam ]]; then
				
				echo "Running samtools quality filtering on $BAM..."

				samtools index -@ 20 $BAM # indexes bam file
				samtools view -@$THREADS -f 0x02 -q 20 -b $BAM > ${BAM%.bam}.filtered.bam # filtering bam file (-q 20 filtering out mapping read quality below 20)
				samtools index -@$THREADS ${BAM%.bam}.filtered.bam # indexes filtered bam file				  
			fi	

			if [[ ! -s $preprocess/$dataset.$sp.$Tissue.$Cond.$rep.bam.bai ]]; then

				echo "Running picard duplicate removal on $BAM..."

				java -jar $PICARD_HOME/picard.jar MarkDuplicates --REMOVE_DUPLICATES true --ASSUME_SORTED true \
					-M ${BAM%.bam}.filtered_dup_metrics.txt \
					-I ${BAM%.bam}.filtered.bam \
					-O $preprocess/$OUTPREX.filtered.DupRem.bam # removes duplicates from bam files

				###### appending replicate name to dataset #######

				if [[ -f $preprocess/$OUTPREX.filtered.DupRem.bam ]]; then

					echo "Appending replicate name to $dataset $sp $Tissue $Cond..."

					mv $preprocess/$OUTPREX.filtered.DupRem.bam $preprocess/$dataset.$sp.$Tissue.$Cond.$rep.bam
					samtools index -@$THREADS $preprocess/$dataset.$sp.$Tissue.$Cond.$rep.bam
                               	fi
			fi

		done
	done
fi

#################################### Quality checking replicate similarity by Species, Tissue, and Condition tags ########################################

if [[ $RUNCHECK == 1 ]]; then

	for dataset in D3 OM; do	
			
		Libs=$(grep -w $sp $in_path/${dataset}SeqData.tsv | cut -f2 | sort -V ) # finds unique library identifiers associated with each species

		if [[ -z $Libs ]]; then
			echo "No libraries found for species $sp"
			continue
		fi

		for lib in $Libs; do # first sections iterates over each library and extracts the tissue, condition, and srrid for each library
			Tissue=$(grep -w $lib $in_path/${dataset}SeqData.tsv | cut -f3 )
			Cond=$(grep -w $lib $in_path/${dataset}SeqData.tsv | cut -f4 )
			rep=$(awk -F'\t' -v lib=$lib '$2==lib {print $6}' $in_path/${dataset}SeqData.tsv | head -n1)	

			cons=$out_path/conservation
			preprocess=$cons/preprocessing
			qc=$preprocess/replicateqc

			mkdir -p $qc
			
			if [[ ! -s $preprocess/$dataset.$sp.$Tissue.$Cond.$rep.coverage.bw ]]; then

				echo "Running bamCoverage on $dataset $sp $Tissue $Cond $rep..."

				bamCoverage --bam $preprocess/$dataset.$sp.$Tissue.$Cond.$rep.bam --normalizeUsing CPM -o $preprocess/$dataset.$sp.$Tissue.$Cond.$rep.coverage.bw
			fi
	
			########################################## quality control for replicate similarity ##########################################
			
			RUNQC=0

			if [[ $RUNQC == 1 ]]; then

				bw_files=($preprocess/*.$sp.*.*.rep*.coverage.bw)
				sorted_bws=$(printf "%s\n" "${bw_files[@]}" | sort -t'.' -k4,4 -k3,3 -k1,1 -k5,5)

				if [[ ! -s $qc/$sp.replicateQC.npz ]]; then
			   
				       echo "Running multiBigwigSummary on $sp..."	

				       multiBigwigSummary bins -b $sorted_bws --binSize 10000 -o $qc/$sp.replicateQC.npz --outRawCounts $qc/$sp.replicateQC.tab	
			       fi

				if [[ ! -s $qc/$sp.replicateQC.heatmap.pdf ]]; then
					
					echo "Plotting heatmap for $sp..."
					
					plotCorrelation --in $qc/$sp.replicateQC.npz --corMethod spearman --whatToPlot heatmap --skipZeros --plotNumbers --colorMap RdYlBu -o $qc/$sp.replicateQC.heatmap.pdf
				fi
			fi
		done
	done
fi

######################################## Running bamPEFragmentSize to assess ATAC-seq fragment read lengths ####################################

if [[ $RUNFRAGMENT == 1 ]]; then

        for dataset in D3 OM; do

		ASSEMBLY=`ls -1 $assembly_path/$sp.assembly.*.fasta`
		number=$(cat $ASSEMBLY | grep -v ">"  | tr -d '\n' | wc -c)

		echo $number

		Libs=$(grep -w $sp $in_path/${dataset}SeqData.tsv | cut -f2 | sort -V)

			if [[ -z $Libs ]]; then
				echo "No libraries found for $sp"
			fi

		for lib in $Libs; do
			Tissue=$(grep -w $lib $in_path/${dataset}SeqData.tsv | cut -f3 )
			Cond=$(grep -w $lib $in_path/${dataset}SeqData.tsv | cut -f4 )
			rep=$(awk -F'\t' -v lib=$lib '$2==lib {print $6}' $in_path/${dataset}SeqData.tsv | head -n1)

                        cons=$out_path/conservation
                        preprocess=$cons/preprocessing
                        qc=$preprocess/noN_atac_fragments

			mkdir -p $qc

			if [[ $rep == rep1 || $rep == rep2 ]]; then

				if [[ ! -s $preprocess/$dataset.$sp.$Tissue.ATAC.$rep.bam.bai ]]; then

					samtools index -@$THREADS $preprocess/$dataset.$sp.$Tissue.ATAC.$rep.bam
				fi

				if [[ ! -s $qc/$dataset.$sp.$Tissue.ATAC.$rep.fragment_size.png ]]; then

					bamPEFragmentSize \
						-b $preprocess/$dataset.$sp.$Tissue.ATAC.$rep.bam \
						-o $qc/$dataset.$sp.$Tissue.ATAC.$rep.fragment_size.png \
						--plotTitle "$dataset $sp $Tissue $rep: ATAC fragment size" \
						--samplesLabel " " --plotFileFormat png --maxFragmentLength 1000 \
						--table $qc/$dataset.$sp.$Tissue.ATAC.$rep.fragment_size_table.tsv \
						--outRawFragmentLengths $qc/$dataset.$sp.$Tissue.ATAC.$rep.raw_fragment_size.tsv
				fi
			fi
		done
        done
fi

########################################################## Peak calling with MACS3 #############################################################

# uses MACS3 for peak calling and then a custom python script to select best narrowPeaks alongside calculation of FRiP scores

if [[ $RUNPEAKCALLING == 1 ]]; then                                                             

	. $APPS2/MACS/venv/bin/activate

	for dataset in D3 OM; do
		
		ASSEMBLY=`ls -1 $assembly_path/$sp.assembly.*.fasta`
		number=$(cat $ASSEMBLY | grep -v ">"  | tr -d '\n' | wc -c)

		echo $number
		
		Libs=$(grep -w $sp $in_path/${dataset}SeqData.tsv | cut -f2 | sort -V ) # finds unique library identifiers associated with each species

			if [[ -z $Libs ]]; then
				echo "No libraries found for species $sp"
				continue
			fi

		for lib in $Libs; do # first sections iterates over each library and extracts the tissue, condition, and srrid for each library
			Tissue=$(grep -w $lib $in_path/${dataset}SeqData.tsv | cut -f3 )
			Cond=$(grep -w $lib $in_path/${dataset}SeqData.tsv | cut -f4 )
			rep=$(awk -F'\t' -v lib=$lib '$2==lib {print $6}' $in_path/${dataset}SeqData.tsv | head -n1)					

			macs3=$out_path/conservation/noN_macs3Data
			preprocess=$out_path/conservation/preprocessing
			qc=$macs3/frip

			mkdir -p $macs3
			mkdir -p $qc

			########################################## ATAC peak calling ##########################################

			if [[ $Cond == ATAC ]]; then

				if [[ ! -s $macs3/$dataset.$sp.$Tissue.$Cond.$rep.narrow_peaks.narrowPeak ]]; then

					echo -e "Running MACS3 ATAC peak calling on $dataset $sp $Tissue $rep..."

					macs3 callpeak \
						-t $preprocess/$dataset.$sp.$Tissue.$Cond.$rep.bam \
						-n $dataset.$sp.$Tissue.$Cond.$rep.narrow --outdir $macs3 \
						-f BAMPE -g $number --nomodel --call-summit -B -q 0.05

					if [[ ! -s $macs3/$dataset.$sp.$Tissue.$Cond.$rep.bestMACS3.narrowPeak ]]; then

						$BIOSCRIPT/SelectBestSummitMACS3narrowPeaks.py -b $macs3/$dataset.$sp.$Tissue.$Cond.$rep.narrow_peaks.narrowPeak \
							| bedtools sort -i - > $macs3/$dataset.$sp.$Tissue.$Cond.$rep.bestMACS3.narrowPeak
					fi
				fi

			fi

			########################################## ChIP peak calling with input control ##########################################

			if [[ $Cond == met ]]; then

				if [[ ! -s $macs3/$dataset.$sp.$Tissue.$Cond.$rep.narrow_peaks.narrowPeak ]]; then

					echo -e "Running MACS3 ChIP ($Cond) peak calling on $dataset $sp $Tissue $rep..."

					macs3 callpeak \
						-t $preprocess/$dataset.$sp.$Tissue.$Cond.$rep.bam \
						-c $preprocess/$dataset.$sp.$Tissue.input.control.bam \
						-n $dataset.$sp.$Tissue.$Cond.$rep.narrow --outdir $macs3 \
						-f BAMPE -g $number --call-summit -B -q 0.05

					if [[ ! -s $macs3/$dataset.$sp.$Tissue.$Cond.$rep.bestMACS3.narrowPeak ]]; then

						$BIOSCRIPT/SelectBestSummitMACS3narrowPeaks.py -b $macs3/$dataset.$sp.$Tissue.$Cond.$rep.narrow_peaks.narrowPeak \
							| bedtools sort -i - > $macs3/$dataset.$sp.$Tissue.$Cond.$rep.bestMACS3.narrowPeak
					fi
				fi
			fi

		        if [[ $Cond == ac ]]; then

                                if [[ ! -s $macs3/$dataset.$sp.$Tissue.$Cond.$rep.bestMACS3.broadPeak ]]; then

                                        echo -e "Running MACS3 ChIP ($Cond) peak calling on $dataset $sp $Tissue $rep..."

                                        macs3 callpeak \
                                                -t $preprocess/$dataset.$sp.$Tissue.$Cond.$rep.bam \
                                                -c $preprocess/$dataset.$sp.$Tissue.input.control.bam \
                                                -n $dataset.$sp.$Tissue.$Cond.$rep.broad --outdir $macs3 \
                                                --broad -f BAMPE -g $number -B -q 0.01

                                        mv $macs3/$dataset.$sp.$Tissue.$Cond.$rep.broad_peaks.broadPeak $macs3/$dataset.$sp.$Tissue.$Cond.$rep.bestMACS3.broadPeak
                                fi
                        fi		
		
			#################### Calculate FRiP scores for ATAC and chip-seq peak calls ##################

			RUNFRIP=1

			if [[ $RUNFRIP == 1 ]]; then

				if [[ $Cond == ATAC || $Cond == ac || $Cond == met ]]; then

					if [[ ! -s $qc/$sp.frip_scores.tsv ]]; then
						echo -e "dataset\tspecies\ttissue\tcondition\treplicate\ttotal_reads\treads_in_peaks\tfrip" > $qc/$sp.frip_scores.tsv
					fi

					if [[ $(awk -v d=$dataset -v s=$sp -v t=$Tissue -v c=$Cond -v r=$rep 'NR>1 && $1==d && $2==s && $3==t && $4==c && $5==r' \
						$qc/$sp.frip_scores.tsv | wc -l) -eq 0 ]]; then

						total=$(samtools view -c $preprocess/$dataset.$sp.$Tissue.$Cond.$rep.bam)
							sort -k1,1 -k2,2n \
							$macs3/$dataset.$sp.$Tissue.$Cond.$rep.bestMACS3.*Peak | \
							bedtools merge -i - > $qc/$dataset.$sp.$Tissue.$Cond.$rep.peaks.merged.bed

						reads_in_peaks=$(samtools view -b $preprocess/$dataset.$sp.$Tissue.$Cond.$rep.bam | \
							bedtools intersect -u -a - -b $qc/$dataset.$sp.$Tissue.$Cond.$rep.peaks.merged.bed | samtools view -c
						)

						frip=$(awk -v n=$reads_in_peaks -v d=$total 'BEGIN { if (d>0) printf "%.5f", n/d; else print "0.00000" }')

						echo -e "${dataset}\t${sp}\t${Tissue}\t${Cond}\t${rep}\t${total}\t${reads_in_peaks}\t${frip}" >> $qc/$sp.frip_scores.tsv
					fi
				fi
			fi
			
			#################### Calculate TSS enrichment scores for ATAC peak calls ##################
		
	                RUNTSS=1	

			if [[ $RUNTSS == 1 ]]; then
				
				if [[ ! -s $qc/$sp.ataqv.json ]]; then
					
					    cut -f1 $assembly_path/$sp.chrom.sizes > $assembly_path/$sp.autosomes.txt

					    ataqv --name $dataset.$sp.$Tissue.$Cond.$rep \
					    	    --peak-file $macs3/$dataset.$sp.$Tissue.$Cond.$rep.bestMACS3.narrowPeak \
					    	    --tss-file $assembly_path/AnnotatedGenomes/$sp.annotation.TSS.bed \
					    	    --autosomal-reference-file $assembly_path/$sp.autosomes.txt $sp \
					    	    $preprocess/$dataset.$sp.$Tissue.$Cond.$rep.bam >> $qc/$sp.ataqv.json
					
				fi
			fi
		done
	done
fi

####################################################### Running quality checks with deeptools ###################################################

# creating plotHeatmaps for all data types alongside a TSS start variant obtained from the data via a custom python script

if [[ $RUNQUALITYCHECK == 1 ]]; then

	for dataset in D3 OM; do

		ASSEMBLY=`ls -1 $assembly_path/$sp.assembly.*.fasta`
		number=$(cat $ASSEMBLY | grep -v ">"  | tr -d '\n' | wc -c)

		echo $number

		Libs=$(grep -w $sp $in_path/${dataset}SeqData.tsv | cut -f2 | sort -V)

			if [[ -z $Libs ]]; then
				echo "No libraries found for $sp"
			fi

		for lib in $Libs; do 
                        Tissue=$(grep -w $lib $in_path/${dataset}SeqData.tsv | cut -f3 )
                        Cond=$(grep -w $lib $in_path/${dataset}SeqData.tsv | cut -f4 )
                        rep=$(awk -F'\t' -v lib=$lib '$2==lib {print $6}' $in_path/${dataset}SeqData.tsv | head -n1)

                        macs3=$out_path/conservation/noN_macs3Data
                        preprocess=$out_path/conservation/preprocessing
			qc=$macs3/heatmaps

			mkdir -p $qc		

			AnnGenomes=/user/work/tk19812/HeliconiniiProject/HeliconGenomeAlignmentAnnotation/UPDATEannotations

			$BIOSCRIPT3/get_TSSandGeneBody.py -i $AnnGenomes/$sp.*.annotation.CAT.gff3 -o $assembly_path/AnnotatedGenomes/$sp.annotation.TSS.bed

			#compute matrix 
		
       			if [[ $Cond == ac || $Cond == met || $Cond == ATAC ]]; then	

				if [[ ! -s $qc/$dataset.$sp.$Tissue.$Cond.$rep.computeMatrix.gz ]]; then

					computeMatrix reference-point --referencePoint center -b 1000 -a 1000 \
						-S $preprocess/$dataset.$sp.$Tissue.$Cond.$rep.coverage.bw \
						-R $macs3/$dataset.$sp.$Tissue.$Cond.$rep.bestMACS3.*Peak \
						--numberOfProcessors 20	-out $qc/$dataset.$sp.$Tissue.$Cond.$rep.computeMatrix.gz
				fi

				if [[ ! -s $qc/TSS.$dataset.$sp.$Tissue.$Cond.$rep.computeMatrix.gz ]]; then

					computeMatrix reference-point --referencePoint TSS -b 3000 -a 3000 \
						-S $preprocess/$dataset.$sp.$Tissue.$Cond.$rep.coverage.bw \
						-R $assembly_path/AnnotatedGenomes/$sp.annotation.TSS.bed \
						--numberOfProcessors 20 -out $qc/TSS.$dataset.$sp.$Tissue.$Cond.$rep.computeMatrix.gz
				fi

				RUNPLOT=1

				if [[ $RUNPLOT == 1 ]]; then

					if [[ ! -s $qc/$dataset.$sp.$Tissue.$Cond.$rep.plotHeatmap.svg ]]; then

						plotHeatmap --matrixFile $qc/$dataset.$sp.$Tissue.$Cond.$rep.computeMatrix.gz \
							--outFileName $qc/$dataset.$sp.$Tissue.$Cond.$rep.plotHeatmap.svg \
							--plotType lines --averageType mean --whatToShow 'plot, heatmap and colorbar' --legendLocation upper-right
					fi
				
					if [[ ! -s $qc/TSS.$dataset.$sp.$Tissue.$Cond.$rep.plotHeatmap.svg ]]; then

						plotHeatmap --matrixFile $qc/TSS.$dataset.$sp.$Tissue.$Cond.$rep.computeMatrix.gz \
							--outFileName $qc/TSS.$dataset.$sp.$Tissue.$Cond.$rep.plotHeatmap.svg \
							--plotType lines --averageType mean --whatToShow 'plot, heatmap and colorbar' --legendLocation upper-right
					fi
				fi
			fi	
		done
	done
fi

# End of data preprocessing

#######################################################################################################################################################################################
##################################################################### Characterisation of enhancers and promoters #####################################################################
#######################################################################################################################################################################################

############################################### Creating consensus sequencing peak files and plotting heatmaps ###################################################

if [[ $RUNCONSENSUS == 1 ]]; then

        macs3=$out_path/conservation/noN_macs3Data
        preprocess=$out_path/conservation/preprocessing
        qc=$macs3/heatmaps

        mkdir -p $qc

        #################### Build replicate-level species-wide peak unions and consensus BAM/bigWig files ####################

        for Cond in ATAC ac met; do
                for rep in rep1 rep2; do

                        echo "Processing $sp $Cond $rep"

                        if [[ $Cond == ac ]]; then
                                peak_type=bestMACS3.broadPeak
                        else
                                peak_type=bestMACS3.narrowPeak
                        fi

			#### peak files ####

                        peak_files=($macs3/D3.$sp.FW.$Cond.$rep.$peak_type $macs3/D3.$sp.HW.$Cond.$rep.$peak_type $macs3/OM.$sp.FW.$Cond.$rep.$peak_type $macs3/OM.$sp.HW.$Cond.$rep.$peak_type)

                        universal=$macs3/$sp.$Cond.universal.$rep.bed
                        universal_tmp=$(mktemp)

                        echo "Building $sp $Cond $rep universal peak set..."

                        cat ${peak_files[@]} | sort -k1,1 -k2,2n | bedtools merge -i - > $universal_tmp
                        mv $universal_tmp $universal
                        echo -n "$sp $Cond $rep universal peaks: "
                        wc -l < $universal

			#### bam files ####

			bam_files=($preprocess/D3.$sp.FW.$Cond.$rep.bam $preprocess/D3.$sp.HW.$Cond.$rep.bam $preprocess/OM.$sp.FW.$Cond.$rep.bam $preprocess/OM.$sp.HW.$Cond.$rep.bam)
                        consensus_bam=$preprocess/$sp.$Cond.consensus.$rep.bam
                        
			if [[ ! -s $consensus_bam ]]; then

                                echo "Merging BAMs for $sp $Cond $rep..."

                                samtools merge -@ $THREADS -o - ${bam_files[@]} | samtools sort -@ $THREADS -o $consensus_bam
                        fi

                        if [[ ! -s $consensus_bam.bai ]]; then

                                echo "Indexing $consensus_bam..."

                                samtools index -@ $THREADS $consensus_bam
                        fi

			#### creating consensus bw's #####

                        consensus_bw=$preprocess/$sp.$Cond.consensus.$rep.bw

                        if [[ ! -s $consensus_bw ]]; then

                                echo "Converting $sp $Cond $rep BAM to bigWig..."

 				bamCoverage --bam $consensus_bam --normalizeUsing CPM --numberOfProcessors $THREADS -o $consensus_bw
                        fi
                done

                # merge completed rep1 + rep2 peak unions

                echo "Building $sp $Cond universal merged peak set..."

                merged=$macs3/$sp.$Cond.universal.merged.bed
                merged_tmp=$(mktemp)

                cat $macs3/$sp.$Cond.universal.rep1.bed $macs3/$sp.$Cond.universal.rep2.bed
	        	| sort -k1,1 -k2,2n | bedtools merge -i - > $merged_tmp

                mv $merged_tmp $merged

                echo -n "$sp $Cond universal merged peaks: "
                wc -l < $merged


                # compute consensus heatmap matrix

                RUNCOMPUTE=1

                if [[ $RUNCOMPUTE == 1 ]]; then

                        matrix=$qc/$sp.$Cond.consensus.computeMatrix.gz

                        if [[ ! -s $matrix ]]; then

                                echo "Creating compute matrix for $sp $Cond..."

                                computeMatrix reference-point --referencePoint center -b 3000 -a 3000 -S \
                                        $preprocess/$sp.$Cond.consensus.rep1.bw $preprocess/$sp.$Cond.consensus.rep2.bw \
                                        -R $merged --numberOfProcessors $THREADS --missingDataAsZero --skipZeros --samplesLabel "rep 1" "rep 2" -O $matrix
                        fi
                fi
        done


	#################### Plots ####################
        
	# peak-size violin plot

        if [[ ! -s $macs3/$sp.violin.pdf ]]; then

                echo "Plotting $sp peak sizes..."

                $BIOSCRIPT2/violinPeakSizeplot.R $macs3/$sp.ATAC.universal.merged.bed $macs3/$sp.ac.universal.merged.bed $macs3/$sp.met.universal.merged.bed \
                        $macs3/$sp.violin.pdf
        fi

        # consensus heatmaps/profiles

        RUNPLOT=1

        if [[ $RUNPLOT == 1 ]]; then

                for Cond in ATAC ac met; do

                        matrix=$qc/$sp.$Cond.consensus.computeMatrix.gz

                        if [[ ! -s $matrix ]]; then
                                echo "Skipping $sp $Cond plots: matrix not found"
                                continue
                        fi

                        if [[ $Cond == ATAC ]]; then
                                colour="Reds"
                        elif [[ $Cond == ac ]]; then
                                colour="Greens"
                        else
                                colour="Blues"
                        fi


                        if [[ ! -s $qc/$sp.$Cond.profile_heatmap.svg ]]; then

                                echo "Plotting $Cond profile heatmap for $sp..."

                                plotProfile -m $matrix --outFileName $qc/$sp.$Cond.profile_heatmap.svg \
                                        --plotType heatmap --perGroup --colors $colour --plotHeight 20 --plotWidth 10 --dpi 300 --legendLocation upper-right
                        fi


                        if [[ ! -s $qc/$sp.$Cond.consensus.profile.svg ]]; then

                                plotProfile --matrixFile $matrix --outFileName $qc/$sp.$Cond.consensus.profile.svg \
			     		--perGroup --averageType mean --plotHeight 10 --plotWidth 10 --dpi 300 --legendLocation upper-right
                        fi
                done
        fi
fi

############################################### Creating cis regulatory element files ###################################################

if [[ $RUNINTERSECT == 1 ]]; then

        for dataset in D3 OM; do
		
		ASSEMBLY=`ls -1 $assembly_path/$sp.assembly.*.fasta`
		number=$(cat $ASSEMBLY | grep -v ">"  | tr -d '\n' | wc -c)

		echo $number

		Libs=$(grep -w $sp $in_path/${dataset}SeqData.tsv | cut -f2 | sort -V)

			if [[ -z $Libs ]]; then
				echo "No libraries found for $sp"
			fi

		for lib in $Libs; do
			Tissue=$(grep -w $lib $in_path/${dataset}SeqData.tsv | cut -f3 )
			Cond=$(grep -w $lib $in_path/${dataset}SeqData.tsv | cut -f4 )
	                rep=$(awk -F'\t' -v lib=$lib '$2==lib {print $6}' $in_path/*SeqData.tsv | head -n1)

			macs3=$out_path/conservation/noN_macs3Data
			cisBeds=$out_path/conservation/cisBeds
		        preprocess=$out_path/conservation/preprocessing

			mkdir -p $cisBeds

			################################# Intersecting cis elements for regulatory mapping ############################
			
			if [[ ! -s $assembly_path/$sp.TSS_1kb.bed ]]; then

				bedtools slop -i $assembly_path/$sp.annotation.TSS.bed -g $assembly_path/$sp.chrom.sizes -b 1000 > $assembly_path/$sp.TSS_1kb.bed
			fi

			RUNCIS=1

			if [[ $RUNCIS == 1 ]]; then

				if [[ $rep == rep1 || $rep == rep2 ]]; then

					if [[ ! -s $cisBeds/$dataset.$sp.$Tissue.$rep.promoters.bed ]]; then

						echo -e "Identifying putative promoters in $dataset $sp $Tissue $rep..."

						bedtools intersect -u \
							-a $macs3/$dataset.$sp.$Tissue.ATAC.$rep.bestMACS3.narrowPeak \
							-b $macs3/$dataset.$sp.$Tissue.ac.$rep.bestMACS3.broadPeak \
						| bedtools intersect -u \
							-a - \
							-b $macs3/$dataset.$sp.$Tissue.met.$rep.bestMACS3.narrowPeak \
						| bedtools intersect -u \
							-a - \
							-b "$assembly_path/$sp.TSS_1kb.bed" \
                                                | sort -k1,1 -k2,2n \
                                                        > $cisBeds/$dataset.$sp.$Tissue.$rep.promoters.bed

					fi

					if [[ ! -s $cisBeds/$dataset.$sp.$Tissue.$rep.enhancers.bed ]]; then

						echo -e "Identifying putative enhancers in $dataset $sp $Tissue $rep..."

						bedtools intersect -u \
							-a $macs3/$dataset.$sp.$Tissue.ATAC.$rep.bestMACS3.narrowPeak \
							-b $macs3/$dataset.$sp.$Tissue.ac.$rep.bestMACS3.broadPeak \
						| bedtools intersect -v \
							-a - \
							-b $macs3/$dataset.$sp.$Tissue.met.$rep.bestMACS3.narrowPeak \
						| bedtools intersect -v \
                                                        -a - \
                                                        -b "$assembly_path/$sp.TSS_1kb.bed" \
						| sort -k1,1 -k2,2n \
							> $cisBeds/$dataset.$sp.$Tissue.$rep.enhancers.bed
					fi
				fi
			fi

			################################ Merging universal cis elements across replicates ################################

			RUNBEDS=1

			if [[ $RUNBEDS == 1 ]]; then
				
			        if [[ $Cond == ATAC ]]; then

					for cis in promoters enhancers; do

						rep1=$cisBeds/$dataset.$sp.$Tissue.rep1.$cis.bed
						rep2=$cisBeds/$dataset.$sp.$Tissue.rep2.$cis.bed

						if [[ -s $rep1 && -s $rep2 ]]; then

							# Reproducible CREs between replicates

							if [[ ! -s $cisBeds/$dataset.$sp.$Tissue.$cis.replicate.bed ]]; then

								echo "Identifying putative $cis shared in both $dataset $Tissue replicates..."

								bedtools intersect -u -f 0.5 -r -a $rep1 -b $rep2 | sort -k1,1 -k2,2n > $cisBeds/$dataset.$sp.$Tissue.$cis.replicate.bed
							fi
					
							# Union of both replicate CRE sets

							if [[ ! -s $cisBeds/$dataset.$sp.$Tissue.$cis.merged.bed ]]; then

								echo "Merging putative $cis from both $dataset $Tissue replicates..."

								cat $rep1 $rep2 | sort -k1,1 -k2,2n | bedtools merge -i - > $cisBeds/$dataset.$sp.$Tissue.$cis.merged.bed
							fi
						fi
					done
				fi
			fi
		done
	done

	for cis in promoters enhancers; do

		if [[ ! -s $cisBeds/$sp.$cis.merged.bed ]]; then

			cat $cisBeds/D3.$sp.FW.$cis.merged.bed $cisBeds/D3.$sp.HW.$cis.merged.bed $cisBeds/OM.$sp.FW.$cis.merged.bed $cisBeds/OM.$sp.HW.$cis.merged.bed \
				| sort -k1,1 -k2,2n | bedtools merge -i - > $cisBeds/$sp.$cis.merged.bed
		fi

		if [[ ! -s $cisBeds/$sp.$cis.replicate.bed ]]; then

			cat $cisBeds/D3.$sp.FW.$cis.replicate.bed $cisBeds/D3.$sp.HW.$cis.replicate.bed $cisBeds/OM.$sp.FW.$cis.replicate.bed $cisBeds/OM.$sp.HW.$cis.replicate.bed \
				| sort -k1,1 -k2,2n | bedtools merge -i - > $cisBeds/$sp.$cis.replicate.bed
		fi
	done
fi

############################################### Quality checking cis regulatory element files ###################################################

if [[ $RUNCISCHECK == 1 ]]; then

        for dataset in D3 OM; do

                ASSEMBLY=`ls -1 $assembly_path/$sp.assembly.*.fasta`
                number=$(cat $ASSEMBLY | grep -v ">"  | tr -d '\n' | wc -c)

                echo $number

                Libs=$(grep -w $sp "$in_path/${dataset}SeqData.tsv" | cut -f2 | sort -V)

                        if [[ -z $Libs ]]; then
                                echo "No libraries found for $sp"
                        fi

                for lib in $Libs; do
                        Tissue=$(grep -w $lib $in_path/${dataset}SeqData.tsv | cut -f3 )
                        Cond=$(grep -w $lib $in_path/${dataset}SeqData.tsv | cut -f4 )
                        rep=$(awk -F'\t' -v lib=$lib '$2==lib {print $6}' $in_path/*SeqData.tsv | head -n1)

                        macs3=$out_path/conservation/noN_macs3Data
                        cisBeds=$out_path/conservation/cisBeds
                        preprocess=$out_path/conservation/preprocessing

			################################ Plotting enhancer and promoter signal intensity for ATAC, H3K27ac, and H3K4me3 ################################

			if [[ $Cond == ATAC || $Cond == ac || $Cond == met ]]; then

				if [[ ! -s $preprocess/$sp.$Cond.consensus.merged.bw ]]; then

					bigwigCompare -b1 $preprocess/$sp.$Cond.consensus.rep1.bw -b2 $preprocess/$sp.$Cond.consensus.rep2.bw \
						--operation mean -o $preprocess/$sp.$Cond.consensus.merged.bw --binSize 10
				fi
				
				if [[ -f $preprocess/$sp.$Cond.consensus.merged.bw ]]; then

					if [[ ! -s $cisBeds/$sp.cis_matrix.gz ]]; then

						computeMatrix reference-point --referencePoint center -b 3000 -a 3000 \
							-S $preprocess/$sp.ATAC.consensus.merged.bw $preprocess/$sp.ac.consensus.merged.bw $preprocess/$sp.met.consensus.merged.bw \
							-R $cisBeds/$sp.promoters.merged.bed $cisBeds/$sp.enhancers.merged.bed -o $cisBeds/$sp.cis_matrix.gz
					fi

					if [[ ! -s $cisBeds/$sp.cis_profile.pdf ]]; then

						plotProfile -m $cisBeds/$sp.cis_matrix.gz -out $cisBeds/$sp.cis_profile.pdf \
							--perGroup --colors red green blue --plotHeight 10 --plotWidth 10 --legendLocation upper-right
					fi
				fi
			fi
                        
			################################ Permutation test for creation of cis regulatory elements ################################

                        permutations=$cisBeds/permutations

                        mkdir -p $permutations

                        RUNSHUFFLE=1

                        if [[ $RUNSHUFFLE == 1 ]]; then

                                if [[ $rep == rep1 || $rep == rep2 ]]; then

                                        shuffleDir=$permutations/$dataset.$sp.$Tissue.ATAC.$rep.shuffled

                                        if [[ ! -s $shuffleDir/Rep.1000.bed ]]; then

                                                echo -e "Creating permutation data for $dataset $sp $Tissue $rep ATAC..."

                                                mkdir -p $shuffleDir

                                                cut -f1-3 $macs3/$dataset.$sp.$Tissue.ATAC.$rep.bestMACS3.narrowPeak | grep -v '^chrom' > $permutations/$dataset.$sp.$Tissue.$rep.tmp

                                                maxTries=1000000

                                                for i in $(seq 1 1000); do

                                                        bedtools shuffle -maxTries $maxTries -chrom -noOverlapping -seed $i \
                                                                -i $permutations/$dataset.$sp.$Tissue.$rep.tmp -g $ASSEMBLY.fai | bedtools sort -i - > $shuffleDir/Rep.$i.bed
                                                done

                                                echo -e "Done"
                                        fi
                                fi
                        fi
			
			RUNPERMUTATION=1

                        if [[ $RUNPERMUTATION == 1 ]]; then

                                if [[ $Cond == ATAC && $rep == rep2 ]]; then

                                        for cis in promoters enhancers; do

                                                mkdir -p $permutations/$dataset.$sp.$Tissue.$cis.merged.shuffled

                                                for i in $(seq 1 1000); do

                                                        outfile=$permutations/$dataset.$sp.$Tissue.$cis.merged.shuffled/Rep.$i.bed

                                                        if [[ ! -s $outfile ]]; then

                                                                if [[ $cis == promoters ]]; then

                                                                        cat \
                                                                                <(
                                                                                        bedtools intersect -u \
                                                                                                -a $permutations/$dataset.$sp.$Tissue.ATAC.rep1.shuffled/Rep.$i.bed \
                                                                                                -b $macs3/$dataset.$sp.$Tissue.ac.rep1.bestMACS3.broadPeak \
                                                                                        | bedtools intersect -u \
                                                                                                -a - \
                                                                                                -b $macs3/$dataset.$sp.$Tissue.met.rep1.bestMACS3.narrowPeak \
                                                                                        | bedtools intersect -u \
                                                                                                -a - \
                                                                                                -b $assembly_path/$sp.TSS_1kb.bed
                                                                                ) \
                                                                                <(
                                                                                        bedtools intersect -u \
                                                                                                -a $permutations/$dataset.$sp.$Tissue.ATAC.rep2.shuffled/Rep.$i.bed \
                                                                                                -b $macs3/$dataset.$sp.$Tissue.ac.rep2.bestMACS3.broadPeak \
                                                                                        | bedtools intersect -u \
                                                                                                -a - \
                                                                                                -b $macs3/$dataset.$sp.$Tissue.met.rep2.bestMACS3.narrowPeak \
                                                                                        | bedtools intersect -u \
                                                                                                -a - \
                                                                                                -b $assembly_path/$sp.TSS_1kb.bed
                                                                                ) \
                                                                        | sort -k1,1 -k2,2n | bedtools merge -i - > $outfile

                                                                elif [[ $cis == enhancers ]]; then

                                                                        cat \
                                                                                <(
                                                                                        bedtools intersect -u \
                                                                                                -a $permutations/$dataset.$sp.$Tissue.ATAC.rep1.shuffled/Rep.$i.bed \
                                                                                                -b $macs3/$dataset.$sp.$Tissue.ac.rep1.bestMACS3.broadPeak \
                                                                                        | bedtools intersect -v \
                                                                                                -a - \
                                                                                                -b $macs3/$dataset.$sp.$Tissue.met.rep1.bestMACS3.narrowPeak \
                                                                                        | bedtools intersect -v \
                                                                                                -a - \
                                                                                                -b $assembly_path/$sp.TSS_1kb.bed
                                                                                ) \
                                                                                <(
                                                                                        bedtools intersect -u \
                                                                                                -a $permutations/$dataset.$sp.$Tissue.ATAC.rep2.shuffled/Rep.$i.bed \
                                                                                                -b $macs3/$dataset.$sp.$Tissue.ac.rep2.bestMACS3.broadPeak \
                                                                                        | bedtools intersect -v \
                                                                                                -a - \
                                                                                                -b $macs3/$dataset.$sp.$Tissue.met.rep2.bestMACS3.narrowPeak \
                                                                                        | bedtools intersect -v \
                                                                                                -a - \
                                                                                                -b $assembly_path/$sp.TSS_1kb.bed
                                                                                ) \
                                                                        | sort -k1,1 -k2,2n | bedtools merge -i - > $outfile
                                                                fi
                                                        fi
                                                done
                                        done
                                fi
                        fi
		done
	done

        if [[ $RUNPERMUTATION == 1 ]]; then

                summary=$permutations/$sp.merged.fold_enrichment.tsv

                if [[ ! -s $summary ]]; then

                        echo -e "species\tcis\tunion_type\tobserved\tperm_mean\tfold_enrichment" > $summary

                        for cis in promoters enhancers; do

                                counts_file=$permutations/$sp.$cis.merged.counts.txt
                                : > $counts_file

                                for i in $(seq 1 1000); do

                                        cat $permutations/D3.$sp.FW.$cis.merged.shuffled/Rep.$i.bed $permutations/D3.$sp.HW.$cis.merged.shuffled/Rep.$i.bed \
                                                $permutations/OM.$sp.FW.$cis.merged.shuffled/Rep.$i.bed $permutations/OM.$sp.HW.$cis.merged.shuffled/Rep.$i.bed \
                                        | sort -k1,1 -k2,2n | bedtools merge -i - | wc -l >> $counts_file
                                done

                                observed=$(wc -l < $cisBeds/$sp.$cis.merged.bed)
                                perm_mean=$(awk '{sum+=$1} END {print sum/NR}' $counts_file)
                                fold=$(awk -v obs=$observed -v mean=$perm_mean 'BEGIN {printf "%.4f", obs/mean}')

                                echo -e "$sp\t$cis\tmerged\t$observed\t$perm_mean\t$fold" >> $summary
                        done
                fi
        fi
fi

########################################## Running CisRegMap.py to assess peak regulatory region locations #####################################

if [[ $RUNREGMAP == 1 ]]; then

        CisRegMapPath=/user/work/wd24922/software/CisRegMapEnv

	for cis in promoters enhancers; do

		ASSEMBLY=`ls -1 $assembly_path/$sp.assembly.*.fasta`
		number=$(cat $ASSEMBLY | grep -v ">"  | tr -d '\n' | wc -c)

		echo $number

		macs3=$out_path/conservation/noN_macs3Data
		cisBeds=$out_path/conservation/cisBeds
		regmap=$out_path/conservation/regmap

		mkdir -p $regmap

		for union in merged replicate; do

			################################ creates shuffle counts for permutations testing  #################################

			RUNSHUFFLE=0

			if [[ $RUNSHUFFLE == 1 ]]; then

				if [[ ! -s $regmap/shuffled/$sp.$cis.$union.shuffled.fofn ]]; then

					echo -e "Creating permutation data for $sp $cis $union..."
					mkdir -p $regmap/shuffled/$sp.$cis.$union.shuffled
					cut -f1-3 $cisBeds/$sp.$cis.$union.bed | grep -v '^chrom' > $cisBeds/$sp.tmp

					maxTries=1000000
					overlapFraction=0.95

					for i in $(seq 1 10000); do
						bedtools shuffle -maxTries $maxTries -f $overlapFraction -chrom -i $cisBeds/$sp.tmp -g $ASSEMBLY.fai | bedtools sort -i - \
							> $regmap/shuffled/$sp.$cis.$union.shuffled/Rep.$i.bed
					done

					ls -1 $regmap/shuffled/$sp.$cis.$union.shuffled/Rep.*.bed | sort -V > $regmap/shuffled/$sp.$cis.$union.shuffled.fofn

					echo -e "Done"
				fi
			fi

			########################################## runs cisregmap slurm array ###########################################

			RUNREGMAP=1

			if [[ $RUNREGMAP == 1 ]]; then
						
				if [[ ! -s $regmap/$sp.$cis.$union.Class.tsv ]]; then

					echo "Running cis regulatory mapping for $sp $cis $union..."

					WORKDIR=$regmap/work/$sp.$cis.$union
					mkdir -p $WORKDIR
					cd $WORKDIR

					python3 -u $BIOSCRIPT2/CisRegMap.py -t $assembly_path/AnnotatedGenomes/$sp.annotation.TSS.bed \
						-i $assembly_path/$sp.annotation.introns.bed -e $assembly_path/$sp.annotation.exons.bed \
						-p $cisBeds/$sp.$cis.$union.bed -r $regmap/shuffled/$sp.$cis.$union.shuffled.fofn

					mv stacked_peak_distribution.pdf $regmap/$sp.$cis.$union.stacked_peak_distribution.pdf
					mv AssGene.tsv $regmap/$sp.$cis.$union.AssGene.tsv
					mv Class.tsv $regmap/$sp.$cis.$union.Class.tsv
				else
					echo "Skipping $sp $cis $union (already exists)..."
				fi
			fi
		done
	done
fi

########################################## Determining association of putative CREs with known wing development genes #####################################

if [[ $RUNGENEANNOTATION == 1 ]]; then

	module load apps/gffread/0.12.8
	module load diamond/2.1.7-iral
	module load seqkit/2.8.2-qs27

	ASSEMBLY=`ls -1 $assembly_path/$sp.assembly.*.fasta`
	number=$(cat $ASSEMBLY | grep -v ">"  | tr -d '\n' | wc -c)

	echo $number
	
	preprocess=$out_path/conservation/preprocessing
	macs3=$out_path/conservation/noN_macs3Data
	cisBeds=$out_path/conservation/cisBeds
	genes=$out_path/conservation/genes
	plots=$genes/plots

	mkdir -p $genes
	mkdir -p $plots

	HmelGFF=$(ls -1 $assembly_path/Hmel.v3.2.annotation.CAT.gff3)
    	HmelGenome=$(ls -1 $assembly_path/Hmel.assembly.*.fasta)
    	SpGFF=$(ls -1 $assembly_path/$sp.*.annotation.CAT.gff3)

        ############################ Mapping orthologs for target species from H.melpomene #############################

	if [[ ! -s $genes/Hmel.proteins.fa ]]; then

		gffread $HmelGFF -g $HmelGenome -y $genes/Hmel.proteins.fa
	fi

    	if [[ ! -s $genes/$sp.proteins.fa ]]; then

		gffread $SpGFF -g $ASSEMBLY -y $genes/$sp.proteins.fa
    	fi

	#for fa in $genes/*.proteins.fa; do
	#    	sed '/^>/! s/\./X/g' "$fa" > "${fa%.fa}.clean.fa"
	#done

	RUNDB=0
		
	if [[ $RUNDB == 1 ]]; then

		diamond makedb --in $genes/Hmel.proteins.clean.fa -d $genes/Hmel_db
		diamond makedb --in $genes/$sp.proteins.clean.fa -d $genes/${sp}_db
	fi

	IDFIXED_FASTA=$genes/Hmel.proteins.clean.idfixed.fa

	if [[ ! -s $IDFIXED_FASTA ]]; then

	    	seqkit replace -p '^transcript:Hmel\.' -r '' $genes/Hmel.proteins.clean.fa > $IDFIXED_FASTA
	fi

	if [[ ! -s $genes/$sp.proteins.clean.idfixed.fa ]]; then

	    	seqkit replace -p ^transcript:$sp\. -r '' $genes/$sp.proteins.clean.fa > $genes/$sp.proteins.clean.idfixed.fa
	fi

	RUNCANDFA=0

	if [[ $RUNCANDFA == 1 ]]; then
		
		rm -f $genes/candidate_queries.fa
		cut -f3 $genes/Hmel_orthologs.tsv | tr ',' '\n' > $genes/target_isoforms.txt
		seqkit grep -n -f $genes/target_isoforms.txt $genes/Hmel.proteins.clean.idfixed.fa > $genes/candidate_queries.fa
	fi	


        ############################ Reciprocal blast search for target species from H.melpomene #############################

	if [[ ! -s $genes/Hmel_vs_${sp}.tsv ]]; then

		diamond blastp -d $genes/${sp}_db.dmnd -q $genes/candidate_queries.fa -o $genes/Hmel_vs_$sp.tsv -f 6 qseqid sseqid bitscore evalue -k 10
	fi

	if [[ ! -s $genes/${sp}_reciprocal_queries.fa ]]; then

		sort -k1,1 -k3,3gr $genes/Hmel_vs_${sp}.tsv | awk '!seen[$1]++' > $genes/Hmel_vs_$sp.best.tsv
		cut -f2 $genes/Hmel_vs_${sp}.best.tsv | sed 's/^transcript://g' | sed "s/^${sp}\.//g" | sort -u > $genes/${sp}_reciprocal_ids.txt
		seqkit grep -n -f $genes/${sp}_reciprocal_ids.txt $genes/${sp}.proteins.clean.idfixed.fa > $genes/${sp}_reciprocal_queries.fa
	fi

	if [[ ! -s $genes/${sp}_vs_Hmel.tsv ]]; then

		diamond blastp -q $genes/${sp}_reciprocal_queries.fa -d $genes/Hmel_db.dmnd -o $genes/${sp}_vs_Hmel.tsv -f 6 qseqid sseqid bitscore evalue
	fi

	if [[ ! -s $genes/${sp}_vs_Hmel.best.tsv ]]; then

	    	sort -k1,1 -k3,3gr $genes/${sp}_vs_Hmel.tsv | awk '!seen[$1]++' > $genes/${sp}_vs_Hmel.best.tsv
	fi

	if [[ ! -s $genes/${sp}_rbh.tsv ]]; then

	    	awk '
	    	function norm(x) {
			gsub(/^transcript:/, "", x)
			sub(/^[A-Za-z]+[.]/, "", x)
			sub(/\.[0-9]+$/, "", x)
			return x
	    	}

	    	NR==FNR {
			q = norm($1)
			s = norm($2)
			best[q] = s
			next
	    	}

	    	{
			q2 = norm($1)
			s2 = norm($2)

			if (best[s2] == q2) {
		    		print $0
			}
	    	}
	    	' \
	    	$genes/Hmel_vs_$sp.best.tsv $genes/${sp}_vs_Hmel.best.tsv > $genes/${sp}_rbh.tsv
	fi

        ############################ Determining species promoter and enhancer lap with target genes #############################

	if [[ ! -s $genes/$sp.genes.bed ]]; then

		awk '$3=="gene" {
	    		match($9, /gene:([^;]+)/, a);
		    	if (a[1] != "") {
				gsub(/^Hhim\./, "", a[1]);
				print $1"\t"$4"\t"$5"\t"a[1]"\t.\t"$7
		    	}
		}' $SpGFF > $genes/$sp.genes.bed
	fi
	
	if [[ ! -s $genes/$sp.rbh.genes.bed || ! -s $genes/$sp.rbh.TSS.bed ]]; then
     
		awk '$3=="transcript" {
			match($9, /transcript:([^;]+)/, a);
			if (a[1] != "") {
		    		gsub(/^[A-Za-z]+\./, "", a[1]);
		    		strand=$7
		    		if (strand == "+") {
					start=$4-1; end=$4
		    		} else {
					start=$5-1; end=$5
	    			}
	    			print $1"\t"start"\t"end"\t"a[1]"\t.\t"strand
			}
    		}' $SpGFF | sed -E 's/\.[0-9]+$//' > $genes/$sp.TSS.bed

	    	cut -f1 $genes/${sp}_rbh.tsv | sed -E 's/\.[0-9]+$//' | sort -u > $genes/${sp}_rbh.genes.txt

	    	grep -F -f $genes/${sp}_rbh.genes.txt $genes/$sp.genes.bed > $genes/$sp.rbh.genes.bed
	    	grep -F -f $genes/${sp}_rbh.genes.txt $genes/$sp.TSS.bed > $genes/$sp.rbh.TSS.bed
	fi

	RUNINTERSECT=0

	if [[ $RUNINTERSECT == 1 ]]; then

		for union in merged replicate; do

			if [[ ! -f $genes/${sp}_rbh.promoters.$union.tsv ]]; then

				bedtools window -a $genes/$sp.rbh.TSS.bed -b $cisBeds/$sp.promoters.$union.bed -w 2000 > $genes/${sp}_rbh.promoters.$union.tsv
			fi

			if [[ ! -s $genes/${sp}_rbh.enhancer.distal.$union.tsv ]]; then

				awk '{
					split($1,a,".");
					print a[1]
			    	}' $genes/${sp}_rbh.tsv | sort -u > $genes/${sp}_rbh_gene_ids.txt

			    	awk -v sp="$sp" '{print "ID=gene:" sp "." $0}' $genes/${sp}_rbh_gene_ids.txt > $genes/${sp}_rbh_gene_ids.prefixed.txt

			    	grep -F -f $genes/${sp}_rbh_gene_ids.prefixed.txt $SpGFF | awk '$3=="gene"' \
			    	| awk '{
					match($9, /ID=gene:([^;]+)/, a);
					print $1"\t"$4"\t"$5"\t"a[1]"\t.\t"$7
			    	}' > $genes/${sp}_rbh_genes.bed

			    	bedtools slop -i $genes/${sp}_rbh_genes.bed -g $assembly_path/$sp.chrom.sizes -b 30000 > $genes/${sp}_rbh_genes.30kb.bed

			    	bedtools intersect -wa -wb -a $genes/${sp}_rbh_genes.30kb.bed -b $cisBeds/$sp.enhancers.$union.bed > $genes/${sp}_rbh.enhancer.distal.$union.tsv
			fi
		done

	fi
	
	RUNVALIDATE=1

	if [[ $RUNVALIDATE == 1 ]]; then
		
		##################### Experimentally validated Lewis2019 optix CREs ########################

		. $APPS2/cactus-bin-v2.9.3/venv-cactus-v2.9.3/bin/activate
		HAL=$assembly_path/79way-Nymphalidae.Cactus-v1.0.hal

		######################## Extract genomes from HAL ########################

		for target in Hlat Hhim Hpet; do

			fasta=$genes/$target.HAL.fasta

			if [[ ! -s $fasta ]]; then

				echo "Extracting $target genome from HAL..."

				hal2fasta $HAL $target > $fasta
			fi

			if [[ ! -s $fasta.fai || $fasta -nt $fasta.fai ]]; then

				echo "Indexing $target HAL FASTA..."

				rm -f $fasta.fai
				samtools faidx $fasta
			fi

			if [[ ! -s $fasta.bwt || $fasta -nt $fasta.bwt ]]; then

				echo "Indexing $target HAL FASTA for BWA..."

				rm -f $fasta.amb $fasta.ann $fasta.bwt $fasta.pac $fasta.sa
				bwa index $fasta
			fi
		done

		######################## Map Lewis2019 sgRNAs to HAL Hlat genome ########################

		if [[ ! -s $genes/Lewis2019_optix_sgRNAs.Hlat.HAL.bed ]]; then

			echo "Mapping Lewis2019 sgRNAs to Hlat..."

			bwa aln -n 1 -o 0 $genes/Hlat.HAL.fasta $genes/Lewis2019_optix_sgRNAs.fa > $genes/Lewis2019_optix_sgRNAs.Hlat.HAL.sai

			bwa samse $genes/Hlat.HAL.fasta $genes/Lewis2019_optix_sgRNAs.Hlat.HAL.sai $genes/Lewis2019_optix_sgRNAs.fa > $genes/Lewis2019_optix_sgRNAs.Hlat.HAL.sam

			samtools view -F 4 $genes/Lewis2019_optix_sgRNAs.Hlat.HAL.sam |
			awk 'BEGIN{OFS="\t"} {
				start=$4-1;
				end=start+length($10);
				strand=and($2,16) ? "-" : "+";
				print $3,start,end,$1,0,strand
			}' > $genes/Lewis2019_optix_sgRNAs.Hlat.HAL.bed
		fi

		######################## Collapse sgRNAs into Hlat CRE intervals ########################

		if [[ ! -s $genes/Lewis2019_optix_CREs.Hlat.HAL.bed ]]; then

			awk 'BEGIN{OFS="\t"} {

				if ($4 ~ /^U1_/) {
					cre="U1"
				}
				else if ($4 ~ /^obs132_/) {
					cre="obs132"
				}
				else if ($4 ~ /^LR1_/) {
					cre="LR1"
				}
				else if ($4 ~ /^LR2_/) {
					cre="LR2"
				}
				else if ($4 ~ /^obs214_/) {
					cre="obs214"
				}
				else {
					next
				}

				if (!(cre in min) || $2 < min[cre])
					min[cre]=$2

				if (!(cre in max) || $3 > max[cre])
					max[cre]=$3

				chr[cre]=$1
			}

			END {
				for (cre in min)
					print chr[cre],min[cre],max[cre],cre
			}' \
			$genes/Lewis2019_optix_sgRNAs.Hlat.HAL.bed | sort -k1,1 -k2,2n > $genes/Lewis2019_optix_CREs.Hlat.HAL.bed
		fi

		######################## Liftover Hlat CREs to Hhim and Hpet ########################

		for target in Hhim Hpet; do

			if [[ ! -s $genes/Lewis2019_optix_CREs.$target.HAL.raw.bed ]]; then

				echo "Lifting Hlat optix CREs to $target..."

				halLiftover $HAL Hlat $genes/Lewis2019_optix_CREs.Hlat.HAL.bed $target $genes/Lewis2019_optix_CREs.$target.HAL.raw.bed
			fi
		done

		######################## Collapse liftover fragments ########################

		# Hlat is already a single interval per CRE
		cp $genes/Lewis2019_optix_CREs.Hlat.HAL.bed $genes/Lewis2019_optix_CREs.Hlat.HAL.spans.bed

		for target in Hhim Hpet; do

			if [[ ! -s $genes/Lewis2019_optix_CREs.$target.HAL.spans.bed ]]; then

				awk 'BEGIN{OFS="\t"} {

					cre=$4

					if (!(cre in min) || $2 < min[cre])
						min[cre]=$2

					if (!(cre in max) || $3 > max[cre])
						max[cre]=$3

					chr[cre]=$1
				}

				END {
					for (cre in min)
						print chr[cre],min[cre],max[cre],cre
				}' $genes/Lewis2019_optix_CREs.$target.HAL.raw.bed | sort -k1,1 -k2,2n > $genes/Lewis2019_optix_CREs.$target.HAL.spans.bed
			fi
		done

		######################## Convert HAL coordinates to local assembly coordinates ########################

		for target in Hlat Hhim Hpet; do

			local_fasta=$(echo $assembly_path/$target.assembly.*.fasta)

			if [[ ! -s $genes/Lewis2019_optix_CREs.$target.HAL.fa ]]; then

				echo "Extracting $target optix CRE sequences..."

				bedtools getfasta -fi $genes/$target.HAL.fasta -bed $genes/Lewis2019_optix_CREs.$target.HAL.spans.bed -nameOnly > $genes/Lewis2019_optix_CREs.$target.HAL.fa
			fi

			if [[ ! -s $local_fasta.bwt ]]; then

				bwa index $local_fasta
			fi

			if [[ ! -s $genes/Lewis2019_optix_CREs.$target.bed ]]; then

				echo "Mapping $target optix CREs to local assembly..."

				bwa mem $local_fasta $genes/Lewis2019_optix_CREs.$target.HAL.fa > $genes/Lewis2019_optix_CREs.$target.local.sam

				samtools view -b -F 2308 -q 20 $genes/Lewis2019_optix_CREs.$target.local.sam | bedtools bamtobed -i - | cut -f1-4 | sort -k1,1 -k2,2n \
					> $genes/Lewis2019_optix_CREs.$target.bed
			fi
		done

		######################## Classify Lewis2019 optix enhancers ########################

		echo -e "CRE\tATAC_bp\tenhancer_bp\tpromoter_bp\tclassification" > $genes/$sp.Lewis2019_optix_CRE_classification.tsv

		for cre in U1 obs132 LR1 LR2 obs214; do

			# CRE absent after liftover
			if ! grep -qw "$cre" $genes/Lewis2019_optix_CREs.$sp.bed; then

				echo -e "$cre\tnot_lifted\tnot_lifted\tnot_lifted\tnot_lifted" >> $genes/$sp.Lewis2019_optix_CRE_classification.tsv

				continue
			fi

			ATAC=$(bedtools intersect -wo -a <(awk -v c=$cre '$4==c' $genes/Lewis2019_optix_CREs.$sp.bed) -b $macs3/$sp.ATAC.universal.merged.bed |	awk '{sum+=$NF} END{print sum+0}')
			enhancer=$(bedtools intersect -wo -a <(awk -v c=$cre '$4==c' $genes/Lewis2019_optix_CREs.$sp.bed) -b $cisBeds/$sp.enhancers.merged.bed | awk '{sum+=$NF} END{print sum+0}')
			promoter=$(bedtools intersect -wo -a <(awk -v c=$cre '$4==c' $genes/Lewis2019_optix_CREs.$sp.bed) -b $cisBeds/$sp.promoters.merged.bed | awk '{sum+=$NF} END{print sum+0}')

			if [[ $enhancer -gt 0 && $promoter -gt 0 ]]; then
				class="mixed enhancer/promoter overlap"
			elif [[ $enhancer -gt 0 ]]; then
				class="enhancer recovered"
			elif [[ $promoter -gt 0 ]]; then
				class="promoter-classified"
			elif [[ $ATAC -gt 0 ]]; then
				class="accessible / uncharacterised"
			else
				class="no detected accessible CRE"
			fi

			echo -e "$cre\t$ATAC\t$enhancer\t$promoter\t$class" >> $genes/$sp.Lewis2019_optix_CRE_classification.tsv

		done

		echo "Lewis2019 optix enhancer classification for $sp:"

		column -t -s $'\t' $genes/$sp.Lewis2019_optix_CRE_classification.tsv

		############################ Plotting promoter and enhancer signal intensity around target genes #############################

		for Cond in ATAC ac met; do

			if [[ ! -s $preprocess/$sp.$Cond.consensus.merged.bam ]]; then
				echo "Merging BAMs for $sp $Cond..."

				samtools merge -@ 20 $preprocess/$sp.$Cond.consensus.merged.bam $preprocess/$sp.$Cond.consensus.rep1.bam $preprocess/$sp.$Cond.consensus.rep2.bam
			fi

			if [[ ! -s $preprocess/$sp.$Cond.consensus.merged.bam.bai ]]; then
				echo "Indexing merged BAM for $sp $Cond..."

				samtools index -@ 20 $preprocess/$sp.$Cond.consensus.merged.bam
			fi

			if [[ ! -s $preprocess/$sp.$Cond.consensus.merged.bw ]]; then
				echo "Generating BigWig for $sp $Cond..."

				bamCoverage --bam $preprocess/$sp.$Cond.consensus.merged.bam --normalizeUsing CPM --numberOfProcessors 20 -o $preprocess/$sp.$Cond.consensus.merged.bw
			fi
		done
	fi
	
	RUNPLOT=1

	if [[ $RUNPLOT == 1 ]]; then
		
		echo "Extracting gene set..."

 		grep -E "G17702|G15037|G9901|G3042|G5855|G8494" $genes/${sp}_rbh_genes.bed > $plots/${sp}_cortex_optix_genes.bed

                bedtools slop -i $plots/${sp}_cortex_optix_genes.bed -g $assembly_path/$sp.chrom.sizes -b 250000 > $plots/${sp}_cortex_optix_windows_250kb.bed

                awk 'BEGIN{OFS="\t"} {print $1":"$2"-"$3, $4}' $plots/${sp}_cortex_optix_windows_250kb.bed  > $plots/${sp}_regions.txt

                echo "Building tracks.ini"

                TRACKS_INI=$plots/${sp}_tracks.ini
                {
                        echo "[ATAC]"
                        echo "file = $preprocess/$sp.ATAC.consensus.merged.bw"
                        echo "title = ATAC"
                        echo "height = 4"
                        echo "color = #E28A8A"
                        echo ""

                        echo "[AC]"
                        echo "file = $preprocess/$sp.ac.consensus.merged.bw"
                        echo "title = AC"
                        echo "height = 4"
                        echo "color = #88CC88"
                        echo ""

                        echo "[MET]"
                        echo "file = $preprocess/$sp.met.consensus.merged.bw"
                        echo "title = MET"
                        echo "height = 4"
                        echo "color = #6988CC"
                        echo ""

                        echo "[Gene body]"
                        echo "file = $plots/${sp}_cortex_optix_genes.bed"
                        echo "title = Genes"
                        echo "height = 1.2"
                        echo "color = grey"
                        echo "alpha = 0.3"
                        echo "border_color = black"
                        echo ""

			echo "[Promoter]"
                        echo "file = $cisBeds/$sp.promoters.merged.bed"
                        echo "title = Experimental Promoters"
                        echo "height = 1"
                        echo "color = black"
                        echo "alpha = 0.4"
                        echo "border_color = black"
                        echo "labels = false"
                        echo "display = collapsed"
                        echo ""

                        echo "[Enhancer]"
                        echo "file = $cisBeds/$sp.enhancers.merged.bed"
                        echo "title = Experimental Enhancers"
                        echo "height = 1"
                        echo "color = #F4A261"
                        echo "alpha = 0.5"
                        echo "border_color = #E76F51"
                        echo "labels = false"
                        echo "display = collapsed"
                        echo ""

                        echo "[Predicted_enhancer]"
                        echo "file = $plots/${sp}_Lewis2019_predicted_enhancers.bed"
                        echo "title = Predicted Enhancers"
                        echo "height = 1"
                        echo "color = #2A9D8F"
                        echo "alpha = 0.7"
                        echo "border_color = #1F6F68"
                        echo "labels = false"
                        echo "display = collapsed"
                        echo ""

			echo "[All_Predicted_enhancer]"
                        echo "file = $plots/$sp.predicted_enhancers.merged.bed"
                        echo "title = All predicted Enhancers"
                        echo "height = 1"
                        echo "color = #2A9D8F"
                        echo "alpha = 0.7"
                        echo "border_color = #1F6F68"
                        echo "labels = false"
                        echo "display = collapsed"
                        echo ""

                        echo "[Lewis2019_CREs]"
                        echo "file = $genes/Lewis2019_optix_CREs.$sp.bed"
                        echo "title = Lewis et al. optix CREs"
                        echo "height = 1"
                        echo "color = #9B5DE5"
                        echo "alpha = 0.6"
                        echo "border_color = #6A4C93"
                        echo "labels = true"
                        echo "display = collapsed"
                        echo ""

		} > $TRACKS_INI
		
		conda activate pygenometracks

                while read region gene; do

                        echo "Plotting $gene at $region"

                        pyGenomeTracks --tracks $TRACKS_INI --region $region --outFileName $plots/${sp}_${gene}_signal_250kb.pdf --dpi 300

                done < $plots/${sp}_regions.txt	
	fi
fi

# End of cis element characterisation

#######################################################################################################################################################################################
####################################################################### Characterisation of transcription factors #####################################################################
#######################################################################################################################################################################################

#################################### Running homer motif analysis to identify overrepresented nucleotide sequences ###############################

if [[ $RUNMOTIF == 1 ]]; then

	for cis in promoters enhancers; do

		ASSEMBLY=`ls -1 $assembly_path/$sp.assembly.*.fasta`
		number=$(cat $ASSEMBLY | grep -v ">"  | tr -d '\n' | wc -c)

		echo $number

		cisBeds=$out_path/conservation/cisBeds
		homer=$out_path/conservation/homer

		mkdir -p $homer

		for union in replicate merged; do

			############################ preparing background from chromosome sizes #############################

                        cisMatrices="$sp.promoters.$union.bed $sp.enhancers.$union.bed"

                        for matrix in $cisMatrices; do
	
				RUNBACKGROUND=1

				if [[ $RUNBACKGROUND == 1 ]]; then

					if [[ ! -s $assembly_path/$sp.chrom.sizes ]]; then

						cut -f1,2 $ASSEMBLY.fai > $assembly_path/$sp.chrom.sizes
					fi

					if [[ -f $assembly_path/$sp.chrom.sizes ]]; then

						base=$(basename $matrix)

						if [[ ! -s $homer/${base%.bed}.bg.bed ]]; then

							bedtools shuffle -i $cisBeds/$matrix -b $assembly_path/$sp.chrom.sizes -chrom > $homer/${base%.bed}.bg.bed
						fi
					fi
				fi

				####################################### promoter and enhancer motif discovery #####################################

				RUNHOMER=1

				if [[ $RUNHOMER == 1 ]]; then

					base=$(basename $matrix)

					if [[ ! -s $homer/${base%.bed}.HomerEnrichedMotif/homerResults.html ]]; then

						$HOMER/findMotifsGenome.pl $cisBeds/$matrix $ASSEMBLY $homer/${base%.bed}.HomerEnrichedMotif \
							-bg $homer/${base%.bed}.bg.bed -mset insects -size given
					fi
				fi
			done
		done
	done	
fi

#################################### Running bedtools nuc to assess gc content of cis-regulatory beds compared to genome ###############################

if [[ $RUNGC == 1 ]]; then

	ASSEMBLY=`ls -1 $assembly_path/$sp.assembly.*.fasta`
	number=$(cat $ASSEMBLY | grep -v ">"  | tr -d '\n' | wc -c)

	echo $number

	cisBeds=$out_path/conservation/cisBeds
	homer=$out_path/conservation/homer

	mkdir -p $homer

	for union in replicate merged; do

		############################ calculating gc content #############################

		cisMatrices="$sp.promoters.$union.bed $sp.enhancers.$union.bed"

		for matrix in $cisMatrices; do

			if [[ ! -s $cisBeds/${matrix%.bed}.GC.txt ]]; then

				echo -e "Calculating GC content for $matrix..."

				bedtools nuc -fi $ASSEMBLY -bed $cisBeds/$matrix > $cisBeds/${matrix%.bed}.GC.txt

				awk '
				NR>1 {
					gc_bases += $5 * $12
					total_bases += $12
				}
				END {
					if (total_bases > 0)
					print "Length-weighted mean GC:", gc_bases/total_bases
				else
					print "No bases found"
				}' $cisBeds/${matrix%.bed}.GC.txt > $cisBeds/${matrix%.bed}.meanGC.txt
			fi

			if [[ ! -s $cisBeds/$sp.assembly.GC.txt ]]; then

				bedtools makewindows -g $assembly_path/$sp.chrom.sizes -w 1000 > $cisBeds/$sp.genome.1kb.windows.bed

				echo -e "Calculating GC content for $ASSEMBLY..."

				bedtools nuc -fi $ASSEMBLY -bed $cisBeds/$sp.genome.1kb.windows.bed > $cisBeds/$sp.assembly.GC.txt

				awk '
				NR>1 {
					effective = $12 - $9
					if (effective > 0) {
						gc_bases += ($5/100) * effective
						total_bases += effective
					}
				}
				END {
			                if (total_bases > 0)
			                        print "Length-weighted mean GC:", gc_bases/total_bases
			                else
		                        print "No valid bases found"
			        }' $cisBeds/$sp.assembly.GC.txt > $cisBeds/$sp.assembly.meanGC.txt
			fi	
		done

		if [[ ! -s $cisBeds/$sp.GC_distribution.pdf ]]; then

			Rscript $BIOSCRIPT2/violinGCplot.R $sp $cisBeds $cisBeds
		fi
	done
fi


################################################### Running motif clustering #############################################################

if [[ $RUNCLUSTER == 1 ]]; then

	ASSEMBLY=`ls -1 $assembly_path/$sp.assembly.*.fasta`
	number=$(cat $ASSEMBLY | grep -v ">"  | tr -d '\n' | wc -c)

	echo $number

	cisBeds=$out_path/conservation/cisBeds
	homer=$out_path/conservation/homer
	cluster=$homer/cluster

	mkdir -p $cluster

	for union in merged replicate; do

		cisMatrices="$sp.promoters.$union.bed $sp.enhancers.$union.bed"

		for matrix in $cisMatrices; do

                        ########################################### Merging known and de novo cis motifs ##########################################

			folder=$homer/${matrix%.bed}.HomerEnrichedMotif

		        known=$folder/knownResults
                        denovo=$folder/homerMotifs.all.motifs


			if [[ ! -s $folder/unified.motifs ]]; then

				cat $denovo $known/known*.motif >> $folder/unified.motifs
			fi

			######################################## Running RSA matrix clustering of cis motifs ######################################

			RUNRSA=1

			if [[ $RUNRSA == 1 ]]; then

				if [[ ! -s $cluster/$folder.cluster ]]; then

					base=$(basename $folder)

					Rscript $RSA/matrix-clustering.R --M $folder/unified.motifs --C ${matrix%.bed} --F homer -o $cluster/$base/$base.cluster --export_heatmap TRUE
				fi
			fi
		done
	done
fi

################################################### Running tobias footprinting and differential binding analysis #########################################################

if [[ $RUNTOBIAS == 1 ]]; then

        ASSEMBLY=`ls -1 $assembly_path/$sp.assembly.*.fasta`
        number=$(cat $ASSEMBLY | grep -v ">"  | tr -d '\n' | wc -c)

        echo $number

	preprocess=$out_path/conservation/preprocessing
        cisBeds=$out_path/conservation/cisBeds
        homer=$out_path/conservation/homer
	tobias=$out_path/conservation/tobias

	mkdir -p $tobias

	if [[ ! -s $preprocess/$sp.ATAC.AllMerged.bam ]]; then

		echo "Merging the following BAMs for $sp..."

		samtools merge -@ 20 -f  - \
			$preprocess/D3.$sp.HW.ATAC.rep1.bam $preprocess/D3.$sp.HW.ATAC.rep2.bam $preprocess/D3.$sp.FW.ATAC.rep1.bam $preprocess/D3.$sp.FW.ATAC.rep2.bam \
			$preprocess/OM.$sp.HW.ATAC.rep1.bam $preprocess/OM.$sp.HW.ATAC.rep2.bam $preprocess/OM.$sp.FW.ATAC.rep1.bam $preprocess/OM.$sp.FW.ATAC.rep2.bam \
			| samtools sort -@ 20 -o $preprocess/$sp.ATAC.AllMerged.bam

		echo "Indexing merged BAM for $sp..."

		samtools index -@ 20 $preprocess/$sp.ATAC.AllMerged.bam
	fi

	conda activate tobias_env

	for union in replicate merged; do

		cisMatrices="$sp.promoters.$union.bed $sp.enhancers.$union.bed"
	
		for matrix in $cisMatrices; do

			RUNFOOTPRINTING=0

			if [[ $RUNFOOTPRINTING == 1 ]]; then

				if [[ ! -s $tobias/${matrix%.bed}_atacorrect.pdf ]]; then

					TOBIAS ATACorrect --bam $preprocess/$sp.ATAC.AllMerged.bam --genome $ASSEMBLY --peaks $cisBeds/$matrix --prefix ${matrix%.bed} --outdir $tobias
				fi

				if [[ ! -s $tobias/${matrix%.bed}_footprintscore.bw ]]; then

					TOBIAS FootprintScores --signal $tobias/${matrix%.bed}_corrected.bw --regions $cisBeds/$matrix --o $tobias/${matrix%.bed}_footprintscore.bw
				fi
			fi

			RUNAGGREGATE=1

			if [[ $RUNAGGREGATE == 1 ]]; then

				if [[ ! -s $tobias/$sp.cis_comparison.$union.pdf ]]; then
					
					echo -e "Plotting signal aggregate comparison for $sp $union..."

					TOBIAS PlotAggregate --TFBS $cisBeds/$sp.promoters.$union.bed $cisBeds/$sp.enhancers.$union.bed --signals $tobias/$sp.promoters.${union}_corrected.bw \
						$tobias/$sp.enhancers.${union}_corrected.bw --output $tobias/$sp.cis_comparison.$union.pdf --share_y both --plot_boundaries --signal-on-x
				fi
			fi

			RUNBINDETECT=1

			if [[ $RUNBINDETECT == 1 ]]; then
				
				if [[ -s $tobias/$sp.promoters.${union}_footprintscore.bw && -s $tobias/$sp.enhancers.${union}_footprintscore.bw ]]; then

					enhancerBed=$cisBeds/$sp.enhancers.$union.bed
					promoterBed=$cisBeds/$sp.promoters.$union.bed
					
					if [[ ! -s $tobias/$sp.cis.$union.bed ]]; then
		
						cat $promoterBed $enhancerBed | sort -k1,1 -k2,2n | bedtools merge > $tobias/$sp.cis.$union.bed
					fi

					if [[ ! -s $tobias/$sp.cis.$union.meme ]]; then

						cat $homer/$sp.promoters.$union.HomerEnrichedMotif/unified.motifs $homer/$sp.enhancers.$union.HomerEnrichedMotif/unified.motifs \
							> $tobias/$sp.cis.$union.motifs
						
						if [[ -f $tobias/$sp.cis.$union.motifs ]]; then

							$BIOSCRIPT2/homer2meme.py $tobias/$sp.cis.$union.motifs $tobias/$sp.cis.$union.meme
						fi
					fi
					
					if [[ ! -s $tobias/$sp.cis.$union.$p.TF/$sp.cis.${union}_results.txt ]]; then

						TOBIAS BINDetect --motifs $tobias/$sp.cis.$union.meme --signals $tobias/$sp.promoters.${union}_footprintscore.bw $tobias/$sp.enhancers.${union}_footprintscore.bw \
							--cond_names Promoter Enhancer 	--genome $ASSEMBLY --peaks $tobias/$sp.cis.$union.bed --prefix $sp.cis.$union --outdir $tobias/$sp.cis.$union.TF					 
					fi
				fi
			fi
		done
	done
fi

################################################### Splitting motifs by TOBIAS BINDetect significance and running TOMTOM to match TFs #########################################################

if [[ $RUNTFMATCH == 1 ]]; then

        ASSEMBLY=`ls -1 $assembly_path/$sp.assembly.*.fasta`
        number=$(cat $ASSEMBLY | grep -v ">"  | tr -d '\n' | wc -c)

        echo $number

        preprocess=$out_path/conservation/preprocessing
        cisBeds=$out_path/conservation/cisBeds
        homer=$out_path/conservation/homer
        tobias=$out_path/conservation/tobias
	tomtom=$out_path/conservation/tomtom

        mkdir -p $tomtom

	cisMatrices="$sp.promoters.merged.bed $sp.enhancers.merged.bed"

	for matrix in $cisMatrices; do

		RUNSPLIT=1

		if [[ $RUNSPLIT == 1 ]]; then
			
			results_file=$tobias/$sp.cis.merged.TF/$sp.cis.merged_results.txt

			if [[ ! -f $results_file ]]; then
				echo "Skipping, results file not found: $results_file"
				continue
			fi

			echo "Processing $results_file..."

			out_prefix=$tomtom/$sp.merged

			awk -v out_prefix=$out_prefix '
			BEGIN { FS=OFS="\t" }

			NR==1 {
				for(i=1;i<=NF;i++){
					if($i=="motif_id") motif_col=i
					if($i=="Promoter_Enhancer_change") change_col=i
					if($i=="Promoter_Enhancer_highlighted") highlight_col=i
				}
				next
			}
			{
				motif=$(motif_col)
				change=$(change_col)+0
				highlight=$(highlight_col)

				# fallback if motif_id missing
				if(motif=="" || motif=="NA" || motif=="-"){
					motif=$1
				}
			}

			highlight=="True" && change<0 { enhancer[motif]=1; next }
			highlight=="True" && change>0 { promoter[motif]=1; next }
			{ nonsig[motif]=1 }

			END {
				for(i in enhancer) print i > (out_prefix ".enhancer.ids")
				for(i in promoter) print i > (out_prefix ".promoter.ids")
				for(i in nonsig)  print i > (out_prefix ".nonsig.ids")
			}
			' $results_file
			
			echo "Created:"
			echo "  ${out_prefix}.enhancer.ids"
			echo "  ${out_prefix}.promoter.ids"
			echo "  ${out_prefix}.nonsig.ids"
		fi
		
		RUNMEME=1				

		if [[ $RUNMEME == 1 ]]; then

			motifdb=$tobias/$sp.cis.merged.meme

			for id_file in $tomtom/$sp.merged.*.ids; do

				type=$(basename $id_file .ids)
				out_meme=$tomtom/$type.meme

				if [[ -s $id_file ]]; then

					awk 'NR==FNR{keep[$1];next}
					!seen && /^MOTIF/{seen=1}
					!seen{print}
					/^MOTIF/{p=($2 in keep)}
					p' $id_file $motifdb > $out_meme

					echo "Created $out_meme"
				fi
			done    	
		fi

		RUNTOMTOM=1

		if [[ $RUNTOMTOM == 1 ]]; then

			if [[ $RUNTOMTOM == 1 ]]; then

				DB=$MOTIF/JASPAR2024_core_insects_non-redundant_meme.txt

				for subset in promoter enhancer nonsig; do

					if [[ ! -s $tomtom/$sp.merged.${subset}_output/tomtom.tsv ]]; then

						echo "Running Tomtom for $tomtom/$sp.merged.$subset.meme"

						tomtom -no-ssc -oc $tomtom/$sp.merged.${subset}_output -min-overlap 5 -dist pearson -evalue -thresh 0.1 $tomtom/$sp.merged.$subset.meme $DB
					fi
				done
			fi
		fi
	done
fi

# End of transcription factor characterisation

#######################################################################################################################################################################################
####################################################################### spatiotemporal conservation with deseq2 #######################################################################
#######################################################################################################################################################################################

############################################################ Pipeline for deseq2 #############################################################

if [[ $RUNDESEQ2 == 1 ]]; then

        ASSEMBLY=`ls -1 $assembly_path/$sp.assembly.*.fasta`
        number=$(cat $ASSEMBLY | grep -v ">"  | tr -d '\n' | wc -c)

        echo $number

	for dataset in D3 OM; do

		Libs=$(grep -w $sp $in_path/${dataset}SeqData.tsv | cut -f2 | sort -V ) # finds unique library identifiers associated with each species

		for lib in $Libs; do # first sections iterates over each library and extracts the tissue, condition, and srrid for each library
			Tissue=$(grep -w $lib $in_path/${dataset}SeqData.tsv | cut -f3 )
			Cond=$(grep -w $lib $in_path/${dataset}SeqData.tsv | cut -f4 )
			rep=$(awk -F'\t' -v lib=$lib '$2==lib {print $6}' $in_path/${dataset}SeqData.tsv | head -n1)	

			cisBeds=$out_path/conservation/cisBeds
			macs3=$out_path/conservation/noN_macs3Data
			preprocess=$out_path/conservation/preprocessing
			deseq2=$out_path/conservation/deseq2
                        genes=$out_path/conservation/genes

			mkdir -p $deseq2

			for union in merged replicate; do
                                        
				######################################## Creates matrix for each signal type and cis element ######################################

				for cis in promoters enhancers; do

					RUNMATRIX=1

					if [[ $RUNMATRIX == 1 ]]; then

						# converting bed to saf for feature counts input

						if [[ ! -s $deseq2/$sp.$cis.$union.saf ]]; then
							
							echo "Converting $sp $cis $union bed to saf..."

							(echo -e "GeneID\tChr\tStart\tEnd\tStrand"; awk 'BEGIN{OFS="\t"} {print $1":"$2+1"-"$3,$1,$2+1,$3,"."}' \
								$cisBeds/$sp.$cis.$union.bed) > $deseq2/$sp.$cis.$union.saf
						fi
						
						RUNCOUNT=0

						if [[ $RUNCOUNT == 1 ]]; then

							if [[ ! -s $deseq2/$sp.$Cond.$cis.$union.Peak-readCountInPeaks.txt ]]; then

								if [[ $Cond == ATAC || $Cond == ac || $Cond == met ]]; then

									if [[ ! -s $deseq2/$sp.$Cond.basename.list ]]; then

										ls $preprocess/*.$sp.*.$Cond.rep[12].bam | sort > $deseq2/$sp.$Cond.list

										while read bam; do
										    	basename "$bam"
										done < $deseq2/$sp.$Cond.list > $deseq2/$sp.$Cond.basename.list
									fi
								fi
						
								if [[ $Cond == ATAC || $Cond == ac || $Cond == met ]]; then
								
									cd $preprocess

									$APPS/featureCounts -p --countReadPairs -a $deseq2/$sp.promoters.$union.saf -F SAF -T 8 \
										-o $deseq2/$sp.$Cond.promoters.$union.Peak-readCountInPeaks.txt $(cat $deseq2/$sp.$Cond.basename.list) \
										2>&1 | tee $deseq2/$sp.$Cond.promoters.$union.featureCounts.log
								fi

								if [[ $Cond == ATAC || $Cond == ac ]]; then

									cd $preprocess

									$APPS/featureCounts -p --countReadPairs -a $deseq2/$sp.enhancers.$union.saf -F SAF -T 8 \
										-o $deseq2/$sp.$Cond.enhancers.$union.Peak-readCountInPeaks.txt $(cat $deseq2/$sp.$Cond.basename.list) \
										2>&1 | tee $deseq2/$sp.$Cond.enhancers.$union.featureCounts.log	
								fi	
							fi
						fi
						
						RUNTSV=0 

						if [[ $RUNTSV == 1 ]]; then

							if [[ ! -s $deseq2/$sp.$Cond.$cis.$union.counts.tsv ]]; then

								if [[ $Cond == ATAC || $Cond == ac || $Cond == met ]]; then

									grep -v '^#' $deseq2/$sp.$Cond.promoters.$union.Peak-readCountInPeaks.txt | cut -f1,7- | sed '1s/\.bam//g' \
										> $deseq2/$sp.$Cond.promoters.$union.counts.tsv
								fi

                                                                if [[ $Cond == ATAC || $Cond == ac ]]; then

									grep -v '^#' $deseq2/$sp.$Cond.enhancers.$union.Peak-readCountInPeaks.txt | cut -f1,7- | sed '1s/\.bam//g' \
                                                                                > $deseq2/$sp.$Cond.enhancers.$union.counts.tsv
								fi
							fi

							if [[ ! -s $deseq2/$sp.$Cond.$cis.$union.coldata.tsv ]]; then
								
								if [[ $Cond == ATAC || $Cond == ac || $Cond == met ]]; then

									echo -e "sample\tdataset\tspecies\ttissue\tcondition\treplicate" > $deseq2/$sp.$Cond.promoters.$union.coldata.tsv

									while read sample; do

										base=${sample%.bam}

										dataset=$(echo $base | cut -d'.' -f1)
										species=$(echo $base | cut -d'.' -f2)
										tissue=$(echo $base | cut -d'.' -f3)
										condition=$(echo $base | cut -d'.' -f4)
										replicate=$(echo $base | cut -d'.' -f5)

										echo -e "$base\t$dataset\t$species\t$tissue\t$condition\t$replicate"

									done < $deseq2/$sp.$Cond.basename.list >> $deseq2/$sp.$Cond.promoters.$union.coldata.tsv
								fi

								if [[ $Cond == ATAC || $Cond == ac ]]; then

									echo -e "sample\tdataset\tspecies\ttissue\tcondition\treplicate" > $deseq2/$sp.$Cond.enhancers.$union.coldata.tsv

                                                                        while read sample; do

                                                                                base=${sample%.bam}

                                                                                dataset=$(echo $base | cut -d'.' -f1)
                                                                                species=$(echo $base | cut -d'.' -f2)
                                                                                tissue=$(echo $base | cut -d'.' -f3)
                                                                                condition=$(echo $base | cut -d'.' -f4)
                                                                                replicate=$(echo $base | cut -d'.' -f5)

                                                                                echo -e "$base\t$dataset\t$species\t$tissue\t$condition\t$replicate"

                                                                        done < $deseq2/$sp.$Cond.basename.list >> $deseq2/$sp.$Cond.enhancers.$union.coldata.tsv
                                                                fi
							fi
						fi

						RUNDESEQ=1

						if [[ $RUNDESEQ == 1 ]]; then

							conda activate deseq2_env

							results=$deseq2/results/$sp
							mkdir -p $results
						
							if [[ $Cond == ATAC || $Cond == ac || $Cond == met ]]; then
									
								if [[ ! -s $deseq2/results/$sp/$sp.$Cond.promoters.$union.PCA.pdf ]]; then

									counts_file=$deseq2/$sp.$Cond.promoters.$union.counts.tsv
									coldata_file=$deseq2/$sp.$Cond.promoters.$union.coldata.tsv

									if [[ -s $counts_file && -s $coldata_file ]]; then

										echo "Running DESeq2:"
										echo "$counts_file"

										GENES_DIR=$genes Rscript $BIOSCRIPT2/deseq2.R $counts_file $coldata_file $results $sp $Cond promoters $union
									else
										echo "Skipping $sp $Cond promoters $union (missing input)..."
									fi
								fi
							fi

							if [[ $Cond == ATAC || $Cond == ac ]]; then
								
								if [[ ! -s $deseq2/results/$sp/$sp.$Cond.enhancers.$union.PCA.pdf ]]; then

									counts_file=$deseq2/$sp.$Cond.enhancers.$union.counts.tsv
									coldata_file=$deseq2/$sp.$Cond.enhancers.$union.coldata.tsv

									if [[ -s $counts_file && -s $coldata_file ]]; then

										echo "Running DESeq2:"
										echo "$counts_file"

										GENES_DIR=$genes Rscript $BIOSCRIPT2/deseq2.R $counts_file $coldata_file $results $sp $Cond enhancers $union
									else
										echo "Skipping $sp $Cond enhancers $union (missing input)..."
									fi
								fi
							fi
						fi
					fi
				done
			done

#			if [[ ! -s $deseq2/results/$sp/${sp}.merged_DE_summary.tsv ]]; then
#
#			    	echo "Merging DE summary tables for $sp"
#
#			    	cat $deseq2/results/$sp/*.DE_summary.tsv \
#				    	| awk 'NR==1 || $0 !~ /^species\t/' \
#				    	> $deseq2/results/$sp/$sp.merged_DE_summary.tsv
#			fi
		done
	done
fi
