library(tidyverse)

# Input/output
input_dir <- "/scratch/project_2019524/pinksalmon_depth_evenness_binned_approx"
output_dir <- "/scratch/project_2019524/pinksalmon_depth_evenness_plots"


dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

samples <- c("8034", "8035", "8036", "8038")

# Read weighted chromosome statistics
files <- list.files(
  input_dir,
  pattern = "_depth_evenness_weighted_chr\\.txt$",
  full.names = TRUE
)

all_stats <- map_dfr(
  files,
  ~ read_tsv(.x, show_col_types = FALSE)
)

# Keep desired samples and methods
plot_data <- all_stats %>%
  filter(
    Sample %in% samples,
    Method %in% c("lowinput", "lowinput_P100", "normalinput")
  ) %>%
  mutate(
    chr_number = as.numeric(chr_number),
    Method = factor(
      Method,
      levels = c("lowinput", "lowinput_P100", "normalinput")
    )
  )

# Make one plot per sample
for (sample_id in samples) {
  
  df <- plot_data %>%
    filter(Sample == sample_id)
  
  p <- ggplot(
    df,
    aes(
      x = factor(chr_number),
      fill = Method,
      color = Method
    )
  ) +
    
    # IQR bar: Q25 to Q75
    geom_linerange(
      aes(
        ymin = Q25_depth,
        ymax = Q75_depth
      ),
      position = position_dodge(width = 0.75),
      linewidth = 5,
      alpha = 0.75
    ) +
    
    # Median
    geom_point(
      aes(
        y = Median_depth
      ),
      position = position_dodge(width = 0.75),
      shape = 16,
      size = 2.5
    ) +
    
    # Mean
    geom_point(
      aes(
        y = Mean_depth
      ),
      position = position_dodge(width = 0.75),
      shape = 4,
      stroke = 1.2,
      size = 3
    ) +
    
    scale_x_discrete(
      drop = FALSE
    ) +
    
    labs(
      title = paste0("Weighted depth distribution — sample ", sample_id),
      subtitle = "Bars = Q25–Q75; filled circles = median; crosses = mean",
      x = "Chromosome",
      y = "Depth",
      fill = "Method",
      color = "Method"
    ) +
    
    theme_classic() +
    
    theme(
      panel.grid.major = element_line(color = "grey85"),
      panel.grid.minor = element_line(color = "grey93"),
      legend.position = "bottom",
      plot.title = element_text(face = "bold"),
      axis.text.x = element_text(angle = 0)
    )
  
  output_file <- file.path(
    output_dir,
    paste0("depth_IQR_mean_median_", sample_id, ".png")
  )
  
  ggsave(
    output_file,
    p,
    width = 12,
    height = 7,
    dpi = 300
  )
  
  message("Written: ", output_file)
}