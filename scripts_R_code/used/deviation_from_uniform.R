#plots how much the distribution of depth 0 bases in the chromosome ends
#deviates from uniform
library(tidyverse)

SAMPLES <- c("8034", "8035", "8036", "8038")

path_input <- "/scratch/project_2019524/pinksalmon_zero_depth_stats"
path_output <- "/scratch/project_2019524/pinksalmon_zero_depth_end_plots"

dir.create(path_output, recursive = TRUE, showWarnings = FALSE)

for (SAMPLE in SAMPLES) {
  
  input_file <- file.path(
    path_input,
    SAMPLE,
    paste0("zero_depth_summary_", SAMPLE, ".csv")
  )
  
  dat <- read_csv(input_file, show_col_types = FALSE)
  
  plot_dat <- dat %>%
    mutate(
      chr_number = as.numeric(chr_number),
      first_15_percent = as.numeric(first_15_percent),
      last_15_percent = as.numeric(last_15_percent)
    ) %>%
    select(
      chr_number,
      chromosome,
      Method,
      first_15_percent,
      last_15_percent
    ) %>%
    pivot_longer(
      cols = c(first_15_percent, last_15_percent),
      names_to = "Region",
      values_to = "Percent"
    ) %>%
    mutate(
      Region = recode(
        Region,
        first_15_percent = "First 15%",
        last_15_percent = "Last 15%"
      ),
      Deviation = Percent - 15
    )
  
  p <- ggplot(
    plot_dat,
    aes(
      x = chr_number,
      y = Deviation,
      color = Method,
      group = interaction(Method, Region),
      linetype = Region,
      shape = Region
    )
  ) +
    geom_hline(
      yintercept = 0,
      linetype = "dashed",
      linewidth = 0.5,
      color = "black"
    ) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2.5) +
    scale_x_continuous(
      breaks = 1:27,
      labels = 1:27
    ) +
    labs(
      title = paste0(
        "Deviation of non-covered bases from 15% — Sample ",
        SAMPLE
      ),
      x = "Chromosome",
      y = "Deviation from 15% (percentage points)",
      color = "Method",
      linetype = "Region",
      shape = "Region"
    ) +
    theme_classic(base_size = 13) +
    theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold"
      )
    )
  
  output_file <- file.path(
    path_output,
    paste0(
      "zero_depth_end_deviation_",
      SAMPLE,
      ".png"
    )
  )
  
  ggsave(
    output_file,
    p,
    width = 12,
    height = 7,
    dpi = 300
  )
  
  message("Saved: ", output_file)
}

