library(tidyverse)
#change the sample IDs
SAMPLES <- c("8034", "8035", "8036", "8038")

#change to the path of the directory that contains the depth evenness/variability stats
path_input <- "/scratch/project_2019524/depth_evenness_stats"
path_output <- "/scratch/project_2019524/depth_evenness_plots"

dir.create(path_output, recursive = TRUE, showWarnings = FALSE)

# fix the 8167 <-> 8036 etc issue
sample_map <- c(
  "8167" = "8036",
  "8168" = "8038"
)

for (SAMPLE in SAMPLES) {
  
  input_ids <- names(sample_map)[sample_map == SAMPLE]
  search_ids <- c(SAMPLE, input_ids)
  
  files <- unlist(lapply(search_ids, function(id) {
    list.files(
      path_input,
      pattern = paste0("^", id, "_.*_depth_evenness\\.txt$"),
      full.names = TRUE
    )
  }))
  
  files <- unique(files)
  
  if (length(files) == 0) {
    warning("No files found for sample ", SAMPLE)
    next
  }
  
  dat <- map_dfr(files, function(file) {
    
    tmp <- read_tsv(file, show_col_types = FALSE)
    
    # Rename the 816_ samples so that we can map them with other 8036 or 8038 samples
    tmp <- tmp %>%
      mutate(
        Sample = recode(
          as.character(Sample),
          "8167" = "8036",
          "8168" = "8038"
        )
      )
    
    tmp
  })
  
  dat <- dat %>%
    mutate(
      Sample = recode(
        as.character(Sample),
        "8167" = "8036",
        "8168" = "8038"
      )
    ) %>%
    filter(Sample == SAMPLE)
  
  plot_dat <- dat %>%
    filter(Level == "Chromosome") %>%
    mutate(
      chr_number = as.numeric(chr_number),
      SD_depth = as.numeric(SD_depth),
      CV_depth = as.numeric(CV_depth)
    ) %>%
    arrange(Method, chr_number)
  
  plot_long <- plot_dat %>%
    select(
      chr_number,
      Method,
      SD_depth,
      CV_depth
    ) %>%
    pivot_longer(
      cols = c(SD_depth, CV_depth),
      names_to = "Metric",
      values_to = "Value"
    ) %>%
    mutate(
      Metric = recode(
        Metric,
        SD_depth = "Standard deviation",
        CV_depth = "Coefficient of variation"
      )
    )
  
  p <- ggplot(
    plot_long,
    aes(
      x = chr_number,
      y = Value,
      color = Method,
      group = Method
    )
  ) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 2.2) +
    facet_wrap(
      ~ Metric,
      scales = "free_y",
      ncol = 1
    ) +
    scale_x_continuous(
      breaks = 1:27,
      labels = 1:27
    ) +
    labs(
      title = paste0("Depth evenness — Sample ", SAMPLE),
      x = "Chromosome",
      y = NULL,
      color = "Method"
    ) +
    theme_classic(base_size = 13) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      strip.text = element_text(face = "bold"),
      legend.position = "bottom",
      
      # Grey grid
      panel.grid.major = element_line(
        color = "grey80",
        linewidth = 0.4
      ),
      panel.grid.minor = element_line(
        color = "grey90",
        linewidth = 0.25
      )
    )
  
  output_file <- file.path(
    path_output,
    paste0("depth_SD_CV_", SAMPLE, ".png")
  )
  
  ggsave(
    output_file,
    p,
    width = 11,
    height = 9,
    dpi = 300
  )
  
  message("Saved: ", output_file)
}

