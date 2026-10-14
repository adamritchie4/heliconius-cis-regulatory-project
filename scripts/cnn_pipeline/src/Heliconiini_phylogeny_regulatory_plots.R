
suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
  library(ape)
})

tree_file <- "MCMCtree.Priors5.Run1-10.Median.Figtree.nex"
cnee_file <- "CNEE_enhancer_overlap_summary.tsv"
out_base <- "All_Species_Regulatory_Distribution_Phylogeny"

zone_levels <- c("Exon", "Intron", "Promoter", "Proximal", "Distal")
zone_colors <- c(Promoter="#6988CC", Proximal="#FFB066", Distal="#88CC88", Exon="#B56272", Intron="#E28A8A")
zone_labels <- c(Exon="Exon", Intron="Intron", Promoter="TSS-proximal", Proximal="Intermediate", Distal="Distal")

species_levels <- c(
  "Hheu", "Htim", "Hpac", "Hcyd", "Hmel", "Hele", "Hatt", "Hhec",
  "Heth", "Hism", "Hbes", "Hnat", "Hwal", "Hbur", "Hxan", "Hdor",
  "Hhie", "Haoe", "Hsap", "Hhew", "Hcon", "Hant", "Hsar", "Hper",
  "Hcha", "Hhim", "Hlat", "Hdem", "Hpet", "Hher", "Htel", "Eisa",
  "Etal", "Ptel"
)

clade_levels <- c(
  "Melpomene", "Wallacei", "Doris", "Aoede", "Sara/Sapho", "Erato", "Eueides", "Other")

clade_map <- c(
  Hheu="Melpomene", Htim="Melpomene", Hpac="Melpomene", Hcyd="Melpomene",
  Hmel="Melpomene", Hatt="Melpomene", Heth="Melpomene",
  Hism="Melpomene", Hbes="Melpomene", Hnat="Melpomene",
  Hwal="Wallacei", Hbur="Wallacei",
  Hxan="Doris", Hdor="Doris", Hhie="Doris",
  Haoe="Aoede",
  Hsap="Sara/Sapho", Hhew="Sara/Sapho", Hcon="Sara/Sapho",
  Hant="Sara/Sapho", Hsar="Sara/Sapho", Hper="Sara/Sapho",
  Hcha="Sara/Sapho", Hele="Sara/Sapho", Hdem="Sara/Sapho",
  Hhim="Erato", Hlat="Erato", Hpet="Erato", Hher="Erato",
  Htel="Erato", Hhec="Erato",
  Eisa="Eueides", Etal="Eueides", Ptel="Other"
)
clade_colors <- c(
  Melpomene="#4778A5", Wallacei="#8166A5", Doris="#3D9787",
  Aoede="#878C91", `Sara/Sapho`="#C98042", Erato="#B15F81",
  Eueides="#879B46", Other="#4C5965"
)

# tree

tree <- ape::read.nexus(tree_file)
if (inherits(tree, "multiPhylo")) tree <- tree[[1]]
species_levels <- intersect(species_levels, tree$tip.label)
tree <- ape::keep.tip(tree, species_levels)
tree <- ape::ladderize(tree, right=TRUE)
tree <- ape::reorder.phylo(tree, order="cladewise")

n_tips <- ape::Ntip(tree)
edge <- tree$edge
x_node <- ape::node.depth.edgelength(tree)
y_node <- ape::node.height(tree)

tip_info <- tibble(Species=tree$tip.label, y=as.numeric(y_node[seq_len(n_tips)]), x=as.numeric(x_node[seq_len(n_tips)])
) %>% arrange(y)

tree_tips <- tip_info %>% mutate(Clade=factor(unname(clade_map[Species]), levels=clade_levels))

horizontal_edges <- tibble(
  x=x_node[edge[, 1]], xend=x_node[edge[, 2]],
  y=y_node[edge[, 2]]
)

parents <- unique(edge[, 1])
vertical_edges <- tibble(
  x=x_node[parents],
  ymin=vapply(parents, function(p) min(y_node[edge[edge[, 1] == p, 2]]), numeric(1)),
  ymax=vapply(parents, function(p) max(y_node[edge[edge[, 1] == p, 2]]), numeric(1))
)

y_bounds <- range(tip_info$y) + c(-0.5, 0.5)

p_tree <- ggplot() +
  geom_segment(
    data=horizontal_edges, aes(x=x, xend=xend, y=y, yend=y),
    linewidth=0.65, colour="#5E6870", lineend="square"
  ) +
  geom_segment(
    data=vertical_edges, aes(x=x, xend=x, y=ymin, yend=ymax),
    linewidth=0.65, colour="#5E6870", lineend="square"
  ) +
  geom_point(data=tree_tips, aes(x=x, y=y, colour=Clade), size=5) +
  geom_point(
    data=filter(tree_tips, Species %in% c("Hhim", "Hlat")),
    aes(x=x, y=y), shape=21, fill=NA, colour="#202630",
    size=7, stroke=1.25
  ) +
  scale_colour_manual(values=clade_colors, breaks=clade_levels, name="Clade") +
  scale_y_continuous(breaks=NULL, expand=expansion(mult=0)) +
  coord_cartesian(ylim=y_bounds, expand=TRUE) +
  theme_classic(base_size=20) +
  theme(
    axis.line.y=element_blank(), axis.ticks.y=element_blank(),
    axis.text.y=element_blank(),
    axis.text.x=element_text(size=20),
    axis.title.x=element_text(size=22),
    plot.margin=margin(8, 12, 8, 12)
  )

if (ape::is.ultrametric(tree, tol=1e-4)) {
  tree_age <- max(x_node[seq_len(n_tips)])
  ages <- pretty(c(0, tree_age), n=4)
  ages <- ages[ages >= 0 & ages <= tree_age]
  p_tree <- p_tree +
    scale_x_continuous(
      breaks=sort(tree_age - ages),
      labels=function(v) scales::number(tree_age - v, accuracy=0.1),
      expand=expansion(mult=c(0.01, 0.30))
    ) + labs(x="Time before present (Ma)")
} else {
  p_tree <- p_tree +
    scale_x_continuous(expand=expansion(mult=c(0.01, 0.20))) +
    labs(x="Evolutionary distance (source tree units)")
}


# read CisRegMap classifications and CNEE overlap table

class_files <- list.files(pattern="^[A-Za-z]{4}\\.Class\\.tsv$")
class_files <- class_files[sub("\\.Class\\.tsv$", "", class_files) %in% tip_info$Species]

all_data <- purrr::map_dfr(class_files, function(f) {
  df <- readr::read_tsv(f, show_col_types=FALSE, progress=FALSE)
  tibble(
    Species=sub("\\.Class\\.tsv$", "", f),
    PeakZone=stringr::str_trim(as.character(df$PeakZone))
  )
})

counts_df <- all_data %>%
  count(Species, PeakZone, name="Count") %>%
  complete(
    Species=tip_info$Species, PeakZone=zone_levels,
    fill=list(Count=0)
  )

percent_df <- counts_df %>%
  group_by(Species) %>%
  mutate(Total=sum(Count),
         Percentage=if_else(Total > 0, Count / Total * 100, 0)) %>%
  ungroup() %>%
  mutate(PeakZone=factor(PeakZone, levels=zone_levels))

totals_df <- percent_df %>% distinct(Species, Total)

cnee_df <- readr::read_tsv(cnee_file, show_col_types=FALSE) %>%
  transmute(
    Species=as.character(species),
    percent_CNEEs_overlap=readr::parse_number(as.character(percent_CNEEs_overlap))
  ) %>%
  filter(Species %in% tip_info$Species, is.finite(percent_CNEEs_overlap))

add_tree_y <- function(df) {
  df %>% mutate(Species=as.character(Species)) %>%
    left_join(tip_info %>% select(Species, y), by="Species") %>%
    filter(is.finite(y))
}
percent_df <- add_tree_y(percent_df)
totals_df <- add_tree_y(totals_df)
cnee_df <- add_tree_y(cnee_df)

bar_y <- function(show_labels=FALSE) {
  scale_y_continuous(
    breaks=if (show_labels) tip_info$y else NULL,
    labels=if (show_labels) tip_info$Species else waiver(),
    expand=expansion(mult=0)
  )
}

bar_theme <- theme_classic(base_size=20) +
  theme(
    axis.text.x=element_text(size=20),
    axis.title.x=element_text(size=22),
    axis.text.y=element_blank(), axis.ticks.y=element_blank(),
    plot.margin=margin(8, 12, 8, 12)
  )

p_count <- ggplot(totals_df, aes(x=Total, y=y)) +
  geom_col(orientation="y", width=0.78, fill="grey85") +
  geom_text(aes(label=scales::comma(Total)), hjust=-0.12, size=7) +
  bar_y(show_labels=TRUE) +
  scale_x_continuous(
    limits=c(0, max(totals_df$Total) * 1.28),
    labels=scales::comma, expand=expansion(mult=0)
  ) +
  coord_cartesian(ylim=y_bounds, expand=TRUE) +
  labs(x="Total predicted enhancers", y=NULL) +
  bar_theme +
  theme(axis.text.y=element_text(size=22, colour="#303740"))

p_percent <- ggplot(percent_df, aes(x=Percentage, y=y, fill=PeakZone)) +
  geom_col(orientation="y", width=0.78) +
  bar_y() +
  scale_fill_manual(
    values=zone_colors, breaks=zone_levels,
    labels=unname(zone_labels[zone_levels]), drop=FALSE
  ) +
  scale_x_continuous(
    breaks=seq(0, 100, 20), labels=function(v) paste0(v, "%"),
    expand=expansion(mult=0)
  ) +
  coord_cartesian(xlim=c(0, 100), ylim=y_bounds, expand=TRUE) +
  labs(x="Genomic distribution", y=NULL) +
  bar_theme

cnee_max <- max(60, ceiling((max(cnee_df$percent_CNEEs_overlap) + 6) / 10) * 10)
p_cnee <- ggplot(cnee_df, aes(x=percent_CNEEs_overlap, y=y)) +
  geom_col(orientation="y", width=0.78, fill="#7886AD") +
  geom_text(
    aes(label=sprintf("%.1f%%", percent_CNEEs_overlap)),
    hjust=-0.12, size=7
  ) +
  bar_y() +
  scale_x_continuous(
    breaks=pretty(c(0, cnee_max), n=4),
    labels=function(v) paste0(v, "%"), expand=expansion(mult=0)
  ) +
  coord_cartesian(xlim=c(0, cnee_max), ylim=y_bounds, expand=FALSE) +
  labs(x="CNEEs overlapping predicted enhancers", y=NULL) +
  bar_theme

p <- (p_tree | plot_spacer() | p_count | plot_spacer() | p_percent | plot_spacer() | p_cnee) +
  plot_layout(ncol=7, widths=c(2.5, 0.20, 3.0, 0.20, 3.0, 0.20, 3.0), guides="collect") +
  plot_annotation(
    title="Predicted enhancer distribution across Heliconiini",
    subtitle="Training species (Hhim and Hlat) indicated by outlined tree tips",
    theme=theme(
      plot.title=element_text(face="bold", size=28),
      plot.subtitle=element_text(size=18)
    )
  ) & theme(legend.position="bottom")

print(p)
ggsave(paste0(out_base, ".png"), p, width=30, height=18, dpi=300, bg="white")
ggsave(paste0(out_base, ".pdf"), p, width=30, height=18, bg="white")
