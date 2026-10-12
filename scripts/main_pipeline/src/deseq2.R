#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)

if(length(args) < 7){
    stop("Usage: Rscript deseq2.R counts.tsv coldata.tsv outdir species condition region union_type")
}

counts_file <- args[1]
coldata_file <- args[2]
outdir <- args[3]
species <- args[4]
condition <- args[5]
region <- args[6]
union_type <- args[7]

suppressPackageStartupMessages({
    library(DESeq2)
    library(ggplot2)
    library(GenomicRanges)
})

cat("====================================\n")
cat("Running DESeq2\n")
cat("Species:", species, "| Condition:", condition, "| Region:", region, "| Union:", union_type, "\n")
cat("====================================\n")

dir.create(outdir, showWarnings=FALSE, recursive=TRUE)

################################################################################
############################ LOAD COUNTS + METADATA ############################
################################################################################

counts <- read.table(counts_file, header=TRUE, row.names=1, sep="\t", check.names=FALSE)
coldata <- read.table(coldata_file, header=TRUE, row.names=1, sep="\t", check.names=FALSE)

if(!all(colnames(counts) %in% rownames(coldata))) stop("Sample mismatch between counts and metadata")
coldata <- coldata[colnames(counts), , drop=FALSE]

if(!all(c("dataset","tissue") %in% colnames(coldata))) stop("Metadata must contain 'dataset' and 'tissue' columns")

coldata$dataset <- factor(coldata$dataset)
coldata$tissue <- factor(coldata$tissue)

if(!"D3" %in% levels(coldata$dataset)) stop("D3 not found in dataset metadata")
if(!"OM" %in% levels(coldata$dataset)) stop("OM not found in dataset metadata")
if(!"FW" %in% levels(coldata$tissue)) stop("FW not found in tissue metadata")
if(!"HW" %in% levels(coldata$tissue)) stop("HW not found in tissue metadata")

coldata$dataset <- relevel(coldata$dataset, ref="D3")
coldata$tissue <- relevel(coldata$tissue, ref="FW")

cat("Samples:", ncol(counts), "\n")
cat("Peaks before filtering:", nrow(counts), "\n")

################################################################################
################################ DESEQ2 CORE ###################################
################################################################################

dds <- DESeqDataSetFromMatrix(countData=counts, colData=coldata, design=~ dataset + tissue)

# Retain peaks with >=10 reads in at least two samples
keep <- rowSums(counts(dds) >= 10) >= 2
dds <- dds[keep, ]

cat("Peaks after filtering:", nrow(dds), "\n")

if(nrow(dds) == 0) stop("No peaks remain after count filtering")

dds <- DESeq(dds)

cat("DESeq2 coefficients:\n")
print(resultsNames(dds))

# Positive LFC = HW > FW or OM > D3
res_tissue <- results(dds, contrast=c("tissue","HW","FW"))
res_dataset <- results(dds, contrast=c("dataset","OM","D3"))

################################################################################
############################ OPTIX/CORTEX LABELS ###############################
################################################################################

cortex_ids <- character(0)
optix_ids <- character(0)

if(region == "promoters"){

    genes_dir <- Sys.getenv("GENES_DIR")
    if(genes_dir == "") warning("GENES_DIR is not set; optix/cortex highlighting will be skipped")

    rbh_file <- file.path(genes_dir, paste0(species, "_rbh.promoters.merged.tsv"))
    saf_file <- file.path(dirname(counts_file), paste0(species, ".", region, ".", union_type, ".saf"))

    if(genes_dir != "" && file.exists(rbh_file) && file.exists(saf_file)){

        rbh <- read.table(rbh_file, header=FALSE, sep="\t", stringsAsFactors=FALSE)
        if(ncol(rbh) < 4) stop("RBH promoter file must contain at least four columns")
        colnames(rbh)[1:4] <- c("Chr","Start","End","GeneID")

        # BED-style start -> GRanges 1-based coordinates
        rbh$Start <- rbh$Start + 1

        rbh_gr <- GRanges(
            seqnames=rbh$Chr,
            ranges=IRanges(start=pmax(1, rbh$Start - 1000), end=rbh$End + 1000),
            gene_id=rbh$GeneID
        )

        saf <- read.table(saf_file, header=TRUE, sep="\t", stringsAsFactors=FALSE)

        required_saf <- c("GeneID","Chr","Start","End")
        if(!all(required_saf %in% colnames(saf))) stop("SAF file must contain GeneID, Chr, Start and End columns")

        saf_gr <- GRanges(seqnames=saf$Chr, ranges=IRanges(start=saf$Start, end=saf$End), peak_id=saf$GeneID)
        hits <- findOverlaps(saf_gr, rbh_gr)

        if(length(hits) > 0){

            peak_ids <- saf_gr$peak_id[queryHits(hits)]
            gene_ids <- mcols(rbh_gr)$gene_id[subjectHits(hits)]

            gene_class_map <- data.frame(
                gene_id=c(
                    "Hhim15G17702.1","Hhim18G15037.1",
                    "Hel_chr15_2G5855.1","Hel_chr18_1G3042.1",
                    "Scf0000031G8494.1","Scf0000026G9901.1"
                ),
                class=c("cortex","optix","cortex","optix","cortex","optix"),
                stringsAsFactors=FALSE
            )

            mapped_class <- gene_class_map$class[match(gene_ids, gene_class_map$gene_id)]

            cortex_ids <- unique(peak_ids[!is.na(mapped_class) & mapped_class == "cortex"])
            optix_ids <- unique(peak_ids[!is.na(mapped_class) & mapped_class == "optix"])
        }

        cat("Cortex promoter peaks:", if(length(cortex_ids)) paste(cortex_ids, collapse=", ") else "none", "\n")
        cat("Optix promoter peaks:", if(length(optix_ids)) paste(optix_ids, collapse=", ") else "none", "\n")

    } else if(genes_dir != "") {
        warning("RBH or SAF file missing; optix/cortex highlighting skipped")
        cat("RBH:", rbh_file, "\n")
        cat("SAF:", saf_file, "\n")
    }
}

################################################################################
################################ MA PLOTS ######################################
################################################################################

plotMA_gg <- function(res, title){

    df <- as.data.frame(res)
    df$peak <- rownames(df)
    df <- df[!is.na(df$padj) & df$baseMean > 0 & !is.na(df$log2FoldChange), , drop=FALSE]

    df$sig <- "NS"
    df$sig[df$padj < 0.05 & df$log2FoldChange > 0] <- "Up"
    df$sig[df$padj < 0.05 & df$log2FoldChange < 0] <- "Down"

    df$category <- df$sig
    df$category[df$peak %in% cortex_ids] <- "Cortex"
    df$category[df$peak %in% optix_ids] <- "Optix"
    df$category <- factor(df$category, levels=c("Up","Down","NS","Cortex","Optix"))

    df_bg <- df[!(df$category %in% c("Cortex","Optix")), , drop=FALSE]
    df_cortex <- df[df$category == "Cortex", , drop=FALSE]
    df_optix <- df[df$category == "Optix", , drop=FALSE]

    ggplot() +
        geom_point(data=df_bg, aes(x=baseMean + 1, y=log2FoldChange, color=category), alpha=0.5, size=0.6) +
        geom_point(data=df_cortex, aes(x=baseMean + 1, y=log2FoldChange, color=category), alpha=1, size=1.2) +
        geom_point(data=df_optix, aes(x=baseMean + 1, y=log2FoldChange, color=category), alpha=1, size=1.2) +
        geom_hline(yintercept=0, linetype="dashed", color="black", alpha=0.4, linewidth=0.4) +
        scale_color_manual(
            values=c("Up"="#E28A8A","Down"="#6988CC","NS"="grey80","Cortex"="#2ca02c","Optix"="#6a3d9a"),
            drop=FALSE
        ) +
        scale_x_log10(name="Mean normalized counts", breaks=c(5,10,50,100,500,1000,5000,10000)) +
        labs(title=title, y="log2 fold change", color=NULL) +
        theme_classic() +
        theme(
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            axis.line=element_line(color="black", linewidth=0.6),
            legend.title=element_blank()
        )
}

################################################################################
################################ OUTPUT MA #####################################
################################################################################

ma_file <- file.path(outdir, paste(species, condition, region, union_type, "MAplots.pdf", sep="."))

pdf(ma_file, width=7, height=6)
print(plotMA_gg(res_tissue, "Tissue HW vs FW"))
print(plotMA_gg(res_dataset, "Dataset OM vs D3"))
dev.off()

################################################################################
################################### PCA ########################################
################################################################################

vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
pca_file <- file.path(outdir, paste(species, condition, region, union_type, "PCA.pdf", sep="."))

pdf(pca_file, width=7, height=6)
print(plotPCA(vsd, intgroup=c("dataset","tissue")))
dev.off()

################################################################################
############################### SAVE RESULTS ###################################
################################################################################

tissue_file <- file.path(outdir, paste(species, condition, region, union_type, "tissue.tsv", sep="."))
dataset_file <- file.path(outdir, paste(species, condition, region, union_type, "dataset.tsv", sep="."))

write.table(as.data.frame(res_tissue), file=tissue_file, sep="\t", quote=FALSE, col.names=NA)
write.table(as.data.frame(res_dataset), file=dataset_file, sep="\t", quote=FALSE, col.names=NA)

cat("====================================\n")
cat("Analysis complete\n")
cat("MA plots:", ma_file, "\n")
cat("PCA:", pca_file, "\n")
cat("Tissue results:", tissue_file, "\n")
cat("Dataset results:", dataset_file, "\n")
cat("====================================\n")
