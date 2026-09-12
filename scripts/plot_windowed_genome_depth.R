#loops through all samples and libraries 
# to produce genome wide depth plots

library(ggplot2)
library(dplyr)
library(readr)
library(tidyr)
library(patchwork)

#change the paths and ids etc accordingly 
INPUT_SAMPLES <- c(
  "8034",
  "8035",
  "8036",
  "8038",
  "8167",
  "8168"
)
#this fixes the 8167/8168 issue
SAMPLE_LABELS <- c(
  "8034" = "8034",
  "8035" = "8035",
  "8036" = "8036",
  "8038" = "8038",
  "8167" = "8036",
  "8168" = "8038"
)

# Methods
METHODS <- c(
  "lowinput",
  "lowinput_P100",
  "normalinput"
)

# base directory containing all the binned coverage data for all samples and methods
#change accordingly
base_dir <- "/scratch/project_2019524/pinksalmon_alignments"
#reference genome path
ref_path <- "/scratch/project_2019524/pinksalmon_reference/GCF_021184085.1_OgorEven_v1.0_genomic.fna.fai"

# save the plots here
output_path <- "/scratch/project_2019524/pinksalmon_depth_plots"


# get chr info from fai

chromosomes <- paste0(
  "NC_",
  sprintf("%06d", 60173:60199),
  ".1"
)

fai <- read.delim(
  ref_path,
  header = FALSE,
  sep = "\t",
  stringsAsFactors = FALSE
)

colnames(fai) <- c(
  "chromosome",
  "chr_length",
  "offset",
  "line_bases",
  "line_width"
)

chrom_lengths <- fai %>%
  filter(chromosome %in% chromosomes) %>%
  mutate(
    chromosome = factor(
      chromosome,
      levels = chromosomes
    ),
    chr_length = as.numeric(chr_length),
    chr_number = match(
      as.character(chromosome),
      chromosomes
    )
  ) %>%
  arrange(chr_number) %>%
  mutate(
    panel = if_else(
      chr_number <= 14,
      "Chromosomes 1–14",
      "Chromosomes 15–27"
    )
  ) %>%
  group_by(panel) %>%
  arrange(chr_number, .by_group = TRUE) %>%
  mutate(
    cumulative_start = lag(
      cumsum(chr_length),
      default = 0
    ),
    cumulative_end = cumsum(chr_length)
  ) %>%
  ungroup()


# label the chrs 
chr_labels <- chrom_lengths %>%
  mutate(
    midpoint = (
      cumulative_start +
        cumulative_end
    ) / 2,
    label = paste0(
      "Chr ",
      chr_number
    )
  )


# name the plots properly

get_method_label <- function(method) {
  
  case_when(
    method == "lowinput" ~ "Low input",
    method == "lowinput_P100" ~ "Low input P100",
    method == "normalinput" ~ "Normal input",
    TRUE ~ method
  )
}

####
####
# plotting function

make_depth_plot <- function(
    input_sample,
    sample_label,
    method
) {
  
  message(
    "Processing sample ",
    input_sample,
    " (label: ",
    sample_label,
    "), method ",
    method
  )
  
  # get the windowed depth file, change path accordingly
  DEPTH_FILE <- paste0(
    base_dir,
    "/",
    method,
    "/coverage_depth/depth_windows/",
    input_sample,
    ".depth_100000bp.txt"
  )
  
  
  # check
  if (!file.exists(DEPTH_FILE)) {
    
    warning(
      "File does not exist: ",
      DEPTH_FILE
    )
    
    return(NULL)
  }
  
  depth <- read.delim(
    DEPTH_FILE,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE
  )
  
  colnames(depth) <- c(
    "chromosome",
    "window_start",
    "window_end",
    "mean_depth"
  )
  
  
  depth <- depth %>%
    mutate(
      chromosome = as.character(chromosome),
      window_start = as.numeric(window_start),
      window_end = as.numeric(window_end),
      mean_depth = as.numeric(mean_depth)
    ) %>%
    left_join(
      chrom_lengths %>%
        mutate(
          chromosome = as.character(chromosome)
        ) %>%
        select(
          chromosome,
          chr_length,
          chr_number,
          panel,
          cumulative_start,
          cumulative_end
        ),
      by = "chromosome"
    ) %>%
    mutate(
      cumulative_position =
        cumulative_start +
        window_start -
        1
    )
  #split the chrs for better readibility
  
  depth_1_14 <- depth %>%
    filter(
      panel == "Chromosomes 1–14"
    )
  
  labels_1_14 <- chr_labels %>%
    filter(
      panel == "Chromosomes 1–14"
    )
  
  boundaries_1_14 <- chrom_lengths %>%
    filter(
      panel == "Chromosomes 1–14",
      cumulative_start > 0
    )
  
  
  p1 <- ggplot(
    depth_1_14,
    aes(
      x = cumulative_position,
      y = mean_depth
    )
  ) +
    
    geom_line(
      linewidth = 0.5,
      color = "grey30",
      alpha = 0.65
    ) +
    
    geom_vline(
      data = boundaries_1_14,
      aes(
        xintercept = cumulative_start
      ),
      linewidth = 0.4,
      linetype = "dotted",
      color = "lightpink3"
    ) +
    
    scale_x_continuous(
      breaks = labels_1_14$midpoint,
      labels = labels_1_14$label,
      expand = c(0.005, 0)
    ) +
    
    scale_y_continuous(
      breaks = seq(
        0,
        35,
        by = 5
      ),
      expand = c(0, 0)
    ) +
    
    coord_cartesian(
      ylim = c(0, 35)
    ) +
    
    labs(
      title = "Chromosomes 1–14",
      x = "Chromosome",
      y = "Mean depth (×)"
    ) +
    
    theme_classic() +
    
    theme(
      plot.title = element_text(
        face = "bold",
        size = 12
      ),
      panel.grid.major.y = element_line(
        color = "grey80"
      ),
      panel.grid.minor.y = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    )
  
  
  # other half of the chrs
  
  depth_15_27 <- depth %>%
    filter(
      panel == "Chromosomes 15–27"
    )
  
  labels_15_27 <- chr_labels %>%
    filter(
      panel == "Chromosomes 15–27"
    )
  
  boundaries_15_27 <- chrom_lengths %>%
    filter(
      panel == "Chromosomes 15–27",
      cumulative_start > 0
    )
  
  
  p2 <- ggplot(
    depth_15_27,
    aes(
      x = cumulative_position,
      y = mean_depth
    )
  ) +
    
    geom_line(
      linewidth = 0.5,
      color = "grey30",
      alpha = 0.65
    ) +
    
    geom_vline(
      data = boundaries_15_27,
      aes(
        xintercept = cumulative_start
      ),
      linewidth = 0.4,
      linetype = "dotted",
      color = "lightpink3"
    ) +
    
    scale_x_continuous(
      breaks = labels_15_27$midpoint,
      labels = labels_15_27$label,
      expand = c(0.005, 0)
    ) +
    
    scale_y_continuous(
      breaks = seq(
        0,
        35,
        by = 5
      ),
      expand = c(0, 0)
    ) +
    
    coord_cartesian(
      ylim = c(0, 35)
    ) +
    
    labs(
      title = "Chromosomes 15–27",
      x = "Chromosome",
      y = "Mean depth (×)"
    ) +
    
    theme_classic() +
    
    theme(
      plot.title = element_text(
        face = "bold",
        size = 12
      ),
      panel.grid.major.y = element_line(
        color = "grey80"
      ),
      panel.grid.minor.y = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    )
  
  
  method_label <- get_method_label(method)
  
  
  final_plot <- p1 / p2 +
    
    plot_annotation(
      title = paste0(
        "Sample ",
        sample_label,
        " — ",
        method_label
      ),
      
      subtitle =
        "Mean sequencing depth across chromosomes",
      
      theme = theme(
        plot.title = element_text(
          face = "bold",
          size = 20,
          hjust = 0.5
        ),
        
        plot.subtitle = element_text(
          size = 12,
          hjust = 0.5
        )
      )
    )
  
  
  return(final_plot)
}


#loop

for (input_sample in INPUT_SAMPLES) {
  
  sample_label <- SAMPLE_LABELS[[input_sample]]
  
  
  sample_output_dir <- file.path(
    output_path,
    sample_label
  )
  
  dir.create(
    sample_output_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  for (method in METHODS) {
    
    final_plot <- make_depth_plot(
      input_sample = input_sample,
      sample_label = sample_label,
      method = method
    )
    
    
    # Skip if input file did not exist
    if (is.null(final_plot)) {
      next
    }
    
    output_file <- file.path(
      sample_output_dir,
      paste0(
        "genome_wide_depth_",
        sample_label,
        "_",
        method,
        ".png"
      )
    )
    
    ggsave(
      filename = output_file,
      plot = final_plot,
      width = 12,
      height = 10,
      units = "in",
      dpi = 300
    )
    
    
    message(
      "Saved: ",
      output_file
    )
  }
}
