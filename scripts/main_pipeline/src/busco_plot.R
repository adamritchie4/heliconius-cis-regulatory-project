library(tidyverse)

busco <- tribble(
  ~Species, ~Category,       ~Count,
  "Hhim",   "Single-copy",    5486,
  "Hhim",   "Duplicated",       95,
  "Hhim",   "Fragmented",       92,
  "Hhim",   "Missing",          87,
  "Hlat",   "Single-copy",    5304,
  "Hlat",   "Duplicated",       44,
  "Hlat",   "Fragmented",       79,
  "Hlat",   "Missing",         333
) %>%
  group_by(Species) %>%
  mutate(Percentage = Count / sum(Count) * 100) %>%
  ungroup()

busco$Category <- factor(busco$Category, levels = c("Single-copy", "Duplicated", "Fragmented", "Missing"))
busco$Species <- factor(busco$Species, levels = c("Hhim", "Hlat"))
busco_cols <- c("Single-copy" = "grey67",  "Duplicated"  = "grey35",  "Fragmented"  = "grey85",  "Missing"     = "#B56272")

p <- ggplot(busco, aes(x = "", y = Count, fill = Category)) +
  geom_col(width = 1, colour = "white", linewidth = 0.3) +
  coord_polar(theta = "y") +
  facet_wrap(~Species, nrow = 1) +
  scale_fill_manual(values = busco_cols) +
  labs(fill = "BUSCO classification") +
  theme_void() +
  theme(strip.text = element_text(face = "bold.italic", size = 12), legend.position = "bottom")

print(p)

ggsave("BUSCO_completeness_pies.png", p, width = 7, height = 4, dpi = 300)