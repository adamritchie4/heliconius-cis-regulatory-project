library(tidyverse)
library(readxl)

files <- c(
  "Hhim.cis.merged_results.xlsx",
  "Hlat.cis.merged_results.xlsx"
)

p_thresh <- 0.05

read_bindetect <- function(f) {
  read_xlsx(f) %>%
    mutate(
      Promoter_Enhancer_pvalue = as.numeric(Promoter_Enhancer_pvalue),
      Promoter_Enhancer_change = as.numeric(Promoter_Enhancer_change),
      highlighted = toupper(as.character(Promoter_Enhancer_highlighted)) %in% c("TRUE", "T", "1", "YES"),
      plot_pvalue = case_when(
        is.na(Promoter_Enhancer_pvalue) ~ NA_real_,
        Promoter_Enhancer_pvalue <= 0 ~ .Machine$double.xmin,
        TRUE ~ Promoter_Enhancer_pvalue
      ),
      neg_log10_p = -log10(plot_pvalue),
      category = case_when(
        highlighted & Promoter_Enhancer_change > 0 ~ "Promoter",
        highlighted & Promoter_Enhancer_change < 0 ~ "Enhancer",
        TRUE ~ "NS"
      )
    )
}

plot_volcano <- function(f) {
  
  df <- read_bindetect(f)
  
  cols <- c(
    "Promoter" = "#6988CC",
    "Enhancer" = "#88CC88",
    "NS" = "grey85"
  )
  
  p <- ggplot(df, aes(x = Promoter_Enhancer_change, y = neg_log10_p)) +
    geom_point(aes(color = category), size = 2, alpha = 0.9) +
    scale_color_manual(values = cols, breaks = c("Promoter", "Enhancer", "NS")) +
    geom_hline(yintercept = -log10(p_thresh), linetype = "dashed", color = "black", linewidth = 0.4) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black", alpha = 0.3, linewidth = 0.4) +
    labs(
      title = tools::file_path_sans_ext(basename(f)),
      x = "Differential binding score",
      y = expression(-log[10](p))
    ) +
    theme_classic() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.6),
      legend.title = element_blank()
    )
  
  outname <- tools::file_path_sans_ext(basename(f))
  
  ggsave(
    paste0(outname, "_volcano.png"),
    p,
    width = 6,
    height = 5,
    dpi = 300
  )
  
  return(p)
}

plots <- lapply(files, plot_volcano)

summarise_file <- function(f) {
  
  df <- read_bindetect(f)
  
  n_total <- nrow(df)
  n_highlighted <- sum(df$highlighted, na.rm = TRUE)
  n_promoter <- sum(df$category == "Promoter", na.rm = TRUE)
  n_enhancer <- sum(df$category == "Enhancer", na.rm = TRUE)
  n_ns <- sum(df$category == "NS", na.rm = TRUE)
  
  n_ns_promoter <- sum(
    df$category == "NS" & df$Promoter_Enhancer_change > 0,
    na.rm = TRUE
  )
  
  n_ns_enhancer <- sum(
    df$category == "NS" & df$Promoter_Enhancer_change < 0,
    na.rm = TRUE
  )
  
  n_p_sig <- sum(df$Promoter_Enhancer_pvalue < p_thresh, na.rm = TRUE)
  
  tibble(
    Species = str_extract(basename(f), "Hhim|Hlat"),
    Dataset = "Merged",
    n_total = n_total,
    
    n_significant = n_highlighted,
    pct_significant = n_highlighted / n_total * 100,
    
    n_p_lt_0.05 = n_p_sig,
    pct_p_lt_0.05 = n_p_sig / n_total * 100,
    
    n_promoter_shift = n_promoter,
    pct_promoter_total = n_promoter / n_total * 100,
    pct_promoter_significant = ifelse(n_highlighted > 0, n_promoter / n_highlighted * 100, NA_real_),
    
    n_enhancer_shift = n_enhancer,
    pct_enhancer_total = n_enhancer / n_total * 100,
    pct_enhancer_significant = ifelse(n_highlighted > 0, n_enhancer / n_highlighted * 100, NA_real_),
    
    promoter_enhancer_ratio_sig = ifelse(n_enhancer > 0, n_promoter / n_enhancer, NA_real_),
    
    n_ns = n_ns,
    pct_ns = n_ns / n_total * 100,
    
    n_ns_promoter_like = n_ns_promoter,
    pct_ns_promoter_like = ifelse(n_ns > 0, n_ns_promoter / n_ns * 100, NA_real_),
    
    n_ns_enhancer_like = n_ns_enhancer,
    pct_ns_enhancer_like = ifelse(n_ns > 0, n_ns_enhancer / n_ns * 100, NA_real_),
    
    ns_direction_bias = mean(
      df$Promoter_Enhancer_change[df$category == "NS"],
      na.rm = TRUE
    ),
    
    ns_abs_effect = mean(
      abs(df$Promoter_Enhancer_change[df$category == "NS"]),
      na.rm = TRUE
    ),
    
    mean_effect = mean(df$Promoter_Enhancer_change, na.rm = TRUE),
    median_effect = median(df$Promoter_Enhancer_change, na.rm = TRUE)
  )
}

summary_table <- map_dfr(files, summarise_file)

print(summary_table, n = Inf, width = Inf)

write_tsv(
  summary_table,
  "species_volcano_summary_table.tsv"
)