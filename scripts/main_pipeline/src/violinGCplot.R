#!/usr/bin/env Rscript


#### script for plotting promoter / enhancer / assembly GC content from text file provided via CisReg.sh
#### Usage: Rscript violinGCplot.R <species> <cisBeds_dir> <output_dir>

suppressPackageStartupMessages({
        library(ggplot2)
        library(dplyr)
})

args <- commandArgs(trailingOnly=TRUE)

if(length(args) < 3){
        stop("Usage: Rscript violinGCplot.R <species> <cisBeds_dir> <output_dir>")
}

sp <- args[1]
cisBeds <- args[2]
outdir <- args[3]

if(!dir.exists(outdir)){
        dir.create(outdir, recursive=TRUE)
}

files <- list(
        promoter=file.path(cisBeds, paste0(sp, ".promoters.merged.GC.txt")),
        enhancer=file.path(cisBeds, paste0(sp, ".enhancers.merged.GC.txt")),
        assembly=file.path(cisBeds, paste0(sp, ".assembly.GC.txt"))
)

read_gc <- function(file, label){

        if(!file.exists(file)) return(NULL)

        df <- read.table(
                file,
                header=TRUE,
                sep="\t",
                comment.char="",
                check.names=FALSE
        )

        gc_col <- grep("pct_gc$", colnames(df), value=TRUE)

        if(length(gc_col) == 0){
                stop(paste("No pct_gc column found in", file))
        }

        data.frame(
                GC=df[[gc_col]],
                Category=label
        )
}

all_data <- bind_rows(
        read_gc(files$promoter, "Promoter"),
        read_gc(files$enhancer, "Enhancer"),
        read_gc(files$assembly, "Assembly")
) %>%
        filter(!is.na(GC))

all_data$Category <- factor(
        all_data$Category,
        levels=c("Assembly","Enhancer","Promoter")
)

summary_table <- all_data %>%
        group_by(Category) %>%
        summarise(
                Mean_GC=mean(GC),
                Median_GC=median(GC),
                SD_GC=sd(GC),
                Min_GC=min(GC),
                Max_GC=max(GC),
                N=n(),
                .groups="drop"
        )

write.table(
        summary_table,
        file=file.path(outdir, paste0(sp, "_GC_summary.tsv")),
        sep="\t",
        quote=FALSE,
        row.names=FALSE
)

cat_colors <- c(
        "Assembly"="grey70",
        "Enhancer"="#88CC88",
        "Promoter"="#6988CC"
)

p <- ggplot(all_data, aes(x=Category, y=GC, fill=Category)) +
        geom_violin(trim=FALSE, alpha=0.8) +
        geom_boxplot(width=0.1, outlier.size=0.3, fill="white") +
        scale_fill_manual(values=cat_colors) +
        labs(
                title=paste("GC Content Distribution -", sp),
                y="GC proportion",
                x=""
        ) +
        theme_minimal() +
        theme(
                legend.position="none",
                panel.border=element_blank(),
                axis.line=element_blank(),
                panel.grid.major.x=element_blank(),
                panel.grid.minor.x=element_blank(),
                panel.grid.minor.y=element_blank(),
                panel.grid.major.y=element_line(color="grey85", linewidth=0.4)
        )

outfile <- file.path(outdir, paste0(sp, "_GC_violin.png"))

ggsave(outfile, p, width=6, height=5, dpi=300)

print(p)
