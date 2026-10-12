library(ggplot2)
library(dplyr)
library(readr)
library(tools)

file_list <- list.files(
  pattern = "FinalTable\\.csv$", 
  full.names = TRUE
)

type_colors <- c("Promoter"="#6988CC", "Enhancer"="#88CC88", "Uncharacterised"="#E28A8A")

combined_data <- lapply(file_list, function(file) {
  df <- read_csv(file)
  
  name_raw <- file_path_sans_ext(basename(file))
  name_clean <- gsub("\\.FinalTable", "", name_raw)   
  name_clean <- gsub("\\.", " ", name_clean)            
  df$Dataset <- name_clean
  
  return(df)
}) %>%
  bind_rows()

summary_df <- combined_data %>%
  filter(Type %in% c("Uncharacterised", "Promoter", "Enhancer")) %>%
  group_by(Dataset, Type) %>%
  summarise(Count = n(), .groups = "drop") %>%
  group_by(Dataset) %>%
  mutate(Percentage = Count / sum(Count) * 100)

plot <- ggplot(summary_df, aes(x = Dataset, y = Percentage, fill = Type)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = paste0(round(Percentage, 1), "%")),
            position = position_stack(vjust = 0.5), size = 3.5) +
  scale_fill_manual(values = type_colors) +
  labs(x = "Dataset", y = "Percentage") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    axis.text.y = element_text(size = 14),
    axis.title = element_text(size = 18),
    legend.title = element_blank(),
    legend.text = element_text(size = 16),
    plot.title = element_text(size = 20, hjust = 0.5),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    panel.grid = element_blank()
  )

ggsave("Combined_PeakType_Distribution.png", plot = plot, width = 12, height = 6, dpi = 300)