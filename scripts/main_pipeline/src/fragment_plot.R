# Load libraries
library(ggplot2)
library(stringr)

# get files
files <- list.files(pattern = "\\.ATAC\\.rep[0-9]+\\.raw_fragment_size\\.tsv$", full.names = TRUE)

# read + clean
read_frag <- function(path) {
  
  df <- read.table(path, sep = "\t", header = TRUE, comment.char = "#", stringsAsFactors = FALSE)
  names(df)[1:2] <- c("Size", "Occurrences")
  
  df$Size <- as.numeric(df$Size)
  df$Occurrences <- as.numeric(df$Occurrences)
  
  fname <- basename(path)
  
  parts <- str_match(fname,
    "^(D3|OM)\\.(Hhim|Hlat|Hpet)\\.(FW|HW)\\.ATAC\\.(rep[0-9]+)\\.raw_fragment_size\\.tsv$"
  )
  
  df$Dataset <- parts[,2]
  df$Species <- parts[,3]
  df$Tissue  <- parts[,4]
  df$Rep     <- parts[,5]
  
  # construct clean sample name
  df$Sample <- paste(df$Dataset, df$Tissue, df$Rep, sep = "_")
  
  return(df)
}

# combine
fragment_data <- do.call(rbind, lapply(files, read_frag))

# normalise per sample
fragment_data$Occurrences <- fragment_data$Occurrences / ave(fragment_data$Occurrences, interaction(fragment_data$Species, fragment_data$Sample), FUN = sum)

# plot
p <- ggplot(fragment_data, aes(x = Size, y = Occurrences, color = Sample, group = interaction(Species, Sample))) +
  geom_line(linewidth = 0.5) +
  facet_wrap(~Species, nrow = 1, scales = "free_y") +
  labs(x = "Fragment Size (bp)", y = "Normalised Count", color = "Sample") +
  
  theme_minimal(base_size = 16) +
  theme(
    strip.text = element_text(size = 16, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 11),
    legend.position = "right"
  ) +
  
  coord_cartesian(xlim = c(0, 800))

# display
print(p)

# save
ggsave("fragment_size_combined_species.png", plot = p, width = 12, height = 5, dpi = 300)