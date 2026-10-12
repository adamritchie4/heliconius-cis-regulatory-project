library(tidyverse)

# ordering definitions
zone_levels  <- c("Exon", "Intron", "Promoter", "Proximal", "Distal")
group_levels <- c("Enhancer", "Promoter")

zone_colors <- c("Promoter"="#6988CC", "Proximal"="#FFB066", "Distal"="#88CC88", "Exon"="#B56272", "Intron"="#E28A8A")

# input files
files <- list.files(pattern = "\\.(enhancers|promoters)\\.merged\\.Class\\.tsv$")

# parse and combine
all_data <- map_dfr(files, function(f) {
  
  df <- read_tsv(f, col_names = c("PeakName", "Locus", "PeakZone"))
  
  tibble(Species  = str_extract(f, "^H\\w+"), Element  = case_when(
      str_detect(f, "enhancers") ~ "Enhancer",
      str_detect(f, "promoters") ~ "Promoter"
    ),
    Dataset  = case_when(str_detect(f, "merged") ~ "Merged",),
    PeakZone = df$PeakZone
  )
})

# count 
counts_df <- all_data %>% 
  count(Species, Dataset, Element, PeakZone, name = "Count") %>%
  mutate(Element  = factor(Element, levels = group_levels),
    PeakZone = factor(PeakZone, levels = zone_levels)
  ) %>%
  complete(Species, Dataset, Element = group_levels, PeakZone = zone_levels, fill = list(Count = 0))

# percentages
percent_df <- counts_df %>%
  group_by(Species, Dataset, Element) %>%
  mutate(Total = sum(Count), Percentage = ifelse(Total > 0, Count / Total * 100, 0)
  ) %>%
  ungroup()

# combinations
combos <- percent_df %>% distinct(Species, Dataset)

# plotting
for (i in seq_len(nrow(combos))) {
  
  sp <- combos$Species[i]
  ds <- combos$Dataset[i]
  
  df_sp <- percent_df %>%
    filter(Species == sp, Dataset == ds) %>%
    mutate(Element  = factor(Element, levels = group_levels), PeakZone = factor(PeakZone, levels = zone_levels)    )
  
  totals_df <- df_sp %>%
    distinct(Element, Total)
  
  p <- ggplot(df_sp, aes(x = Element, y = Percentage, fill = PeakZone)) +
    
    geom_bar(stat = "identity") +
    scale_x_discrete(drop = FALSE) +
    scale_fill_manual(values = zone_colors, breaks = zone_levels, drop = FALSE, na.translate = FALSE) +
    geom_text(data = totals_df, aes(x = Element, y = 105, label = Total), inherit.aes = FALSE, size = 3.5) +
    labs(
      title = paste("Regulatory Region Distribution (%) —", sp, ds),
      x     = NULL,
      y     = "Percentage of Peaks",
      fill  = "Region"
    ) +
        scale_y_continuous(limits = c(0, 110), breaks = seq(0, 100, by = 20), expand = c(0, 0)) +
    theme_classic() +
    theme(
      plot.title.position = "plot",
      plot.margin         = margin(t = 20, r = 10, b = 10, l = 10),
      axis.text.x         = element_text(angle = 0),
      plot.title          = element_text(face = "bold", size = 14),
      legend.position     = "right"
    )
  
  print(p)
  
  ggsave(filename = paste0(sp, "_", ds, "_Regulatory_Stacked_Percent.png"), plot = p, width = 6, height = 5)
}