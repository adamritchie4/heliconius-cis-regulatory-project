library(tidyverse)
library(patchwork)

# definitions

zone_levels <- c("Exon", "Intron", "Promoter", "Proximal", "Distal")
zone_colors <- c("Promoter"="#6988CC", "Proximal"="#FFB066", "Distal"="#88CC88", "Exon"="#B56272", "Intron"="#E28A8A")

species_levels <- c(
  "Hheu", "Htim", "Hpac", "Hcyd", "Hmel", "Hele", "Hatt", "Hhec", "Heth", "Hnum", "Hism", "Hbes", "Hnat",
  "Hwal", "Hbur",
  "Hxan", "Hdor", "Hhie",
  "Haoe",
  "Hsap", "Hhew", "Hcon", "Hant", "Hsar", "Hper", "Hcha",
  "Hhim", "Hlat", "Hdem", "Hpet", "Hher", "Htel",
  "Eisa", "Etal",
  "Ptel"
)

clade_levels <- c("Melpomene", "Wallacei", "Doris", "Aoede", "Sara/Sapho", "Erato", "Eueides", "Other")

clade_map <- c(
  Hheu="Melpomene", Htim="Melpomene", Hpac="Melpomene", Hcyd="Melpomene", Hmel="Melpomene",
  Hele="Melpomene", Hatt="Melpomene", Hhec="Melpomene", Heth="Melpomene", Hnum="Melpomene",
  Hism="Melpomene", Hbes="Melpomene", Hnat="Melpomene",
  Hwal="Wallacei", Hbur="Wallacei",
  Hxan="Doris", Hdor="Doris", Hhie="Doris",
  Haoe="Aoede",
  Hsap="Sara/Sapho", Hhew="Sara/Sapho", Hcon="Sara/Sapho", Hant="Sara/Sapho",
  Hsar="Sara/Sapho", Hper="Sara/Sapho", Hcha="Sara/Sapho",
  Hhim="Erato", Hlat="Erato", Hdem="Erato", Hpet="Erato", Hher="Erato", Htel="Erato",
  Eisa="Eueides", Etal="Eueides",
  Ptel="Other"
)

# read data

files <- list.files(pattern = "^[A-Za-z]{4}\\.Class\\.tsv$")

all_data <- map_dfr(files, function(f) {
  df <- read_tsv(f, show_col_types = FALSE)
  tibble(Species = sub("\\.Class\\.tsv$", "", f), PeakZone = str_trim(df$PeakZone))
})

# check classifications

unknown <- all_data %>%
  filter(!PeakZone %in% zone_levels) %>%
  count(Species, PeakZone)

if (nrow(unknown) > 0) {
  print(unknown)
  stop("Unrecognised PeakZone categories detected")
}

# counts + percentages

counts_df <- all_data %>%
  count(Species, PeakZone, name = "Count") %>%
  complete(Species, PeakZone = zone_levels, fill = list(Count = 0))

percent_df <- counts_df %>%
  group_by(Species) %>%
  mutate(Total = sum(Count), Percentage = Count / Total * 100) %>%
  ungroup() %>%
  mutate(
    Clade = factor(clade_map[Species], levels = clade_levels),
    Species = factor(Species, levels = rev(species_levels)),
    PeakZone = factor(PeakZone, levels = zone_levels)
  )

totals_df <- percent_df %>%
  distinct(Species, Clade, Total)

cnee_df <- read_tsv("CNEE_enhancer_overlap_summary.tsv", show_col_types = FALSE) %>%
  mutate(
    Clade = factor(clade_map[species], levels = clade_levels),
    Species = factor(species, levels = rev(species_levels))
  )

# QC

print(percent_df %>% group_by(Species) %>% summarise(Percentage = sum(Percentage)))

# total count plot

p_count <- ggplot(totals_df, aes(x = Total, y = Species)) +
  geom_col(width = 0.8, fill = "grey85") +
  geom_text(aes(label = scales::comma(Total)), hjust = -0.1, size = 3) +
  facet_grid(Clade ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_x_continuous(labels = scales::comma, expand = expansion(mult = c(0, 0.22))) +
  labs(x = "Total Predicted Enhancers", y = NULL) +
  theme_classic() +
  theme(
    axis.text.y = element_text(size = 9),
    strip.placement = "outside",
    strip.background = element_blank(),
    strip.text.y.left = element_text(face = "bold", angle = 0, size = 10),
    panel.spacing.y = unit(0.6, "lines"),
    plot.margin = margin(10, 5, 10, 10)
  )

# percentage plot

p_percent <- ggplot(percent_df, aes(x = Percentage, y = Species, fill = PeakZone)) +
  geom_col(width = 0.8) +
  facet_grid(Clade ~ ., scales = "free_y", space = "free_y") +
  scale_fill_manual(values = zone_colors, breaks = zone_levels, drop = FALSE) +
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 20), expand = c(0, 0)) +
  labs(x = "Percentage of Predicted Enhancers", y = NULL, fill = "Region") +
  theme_classic() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    strip.background = element_blank(),
    strip.text.y = element_blank(),
    legend.position = "right",
    panel.spacing.y = unit(0.6, "lines"),
    plot.margin = margin(10, 10, 10, 5)
  )

#  CNEE overlap plot

p_cnee <- ggplot(cnee_df, aes(x = percent_CNEEs_overlap, y = Species)) +
  geom_col(width = 0.8, fill = "#7A9E7E") +
  geom_text(aes(label = paste0(percent_CNEEs_overlap, "%")), hjust = -0.1, size = 3) +
  facet_grid(Clade ~ ., scales = "free_y", space = "free_y") +
  scale_x_continuous(limits = c(0, 60), breaks = seq(0, 60, 20), expand = c(0, 0)) +
  labs(x = "CNEEs Overlapping Predicted Enhancers (%)", y = NULL) +
  theme_classic() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    strip.background = element_blank(),
    strip.text.y = element_blank(),
    panel.spacing.y = unit(0.6, "lines"),
    plot.margin = margin(10, 10, 10, 5)
  )

# combine

p <- p_count + p_percent + p_cnee +
  plot_layout(widths = c(1.1, 2.4, 1.1)) +
  plot_annotation(
    title = "Regulatory Region Distribution",
    theme = theme(plot.title = element_text(face = "bold", size = 14))
  )

print(p)

ggsave("All_Species_Regulatory_Distribution_Counts_CNEE.png",  p,  width = 17,  height = 14,  dpi = 300)