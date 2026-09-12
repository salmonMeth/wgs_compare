## creates GC content and mapping quality score distribution plots
# using the stats.txt files
library(tidyverse)


# change the paths, sample IDs and the method accordingly
BAMDIR <- "/scratch/project_2019524/pinksalmon_alignments/lowinput_P100"
SAMPLES <- c("8034", "8035", "8036", "8038")
METHOD <- "Low Input-P100"

# Output directory
OUTDIR <- file.path(BAMDIR, "plots")

dir.create(
  OUTDIR,
  recursive = TRUE,
  showWarnings = FALSE
)
#get the stats files
stats_files <- list.files(
  path = BAMDIR,
  pattern = "\\.stats\\.txt$",
  recursive = TRUE,
  full.names = TRUE
)

if (length(stats_files) == 0) {
  stop(
    "No .stats.txt files found under: ",
    BAMDIR
  )
}

cat("Found", length(stats_files), "stats files:\n")

print(stats_files)


read_stats <- function(sample) {
  
  file <- stats_files[
    SAMPLES == sample
  ]
  
  if (length(file) == 0) {
    stop(
      "No stats file found for sample: ",
      sample
    )
  }
  
  readLines(file)
}


#extract the gc content 
extract_gc <- function(sample) {
  
  lines <- read_stats(sample)
  
  # GCF = GC content of first fragments / R1
  gcf <- lines[
    grepl("^GCF\\t", lines)
  ]
  
  # GCL = GC content of last fragments / R2
  gcl <- lines[
    grepl("^GCL\\t", lines)
  ]
  
  if (length(gcf) == 0) {
    warning(
      "No GCF section found for sample ",
      sample
    )
  }
  
  if (length(gcl) == 0) {
    warning(
      "No GCL section found for sample ",
      sample
    )
  }
  #read 1
  gcf_df <- NULL
  
  if (length(gcf) > 0) {
    
    gcf_df <- read_tsv(
      I(paste(gcf, collapse = "\n")),
      col_names = c(
        "type",
        "gc",
        "count"
      ),
      show_col_types = FALSE
    ) %>%
      mutate(
        sample = sample,
        read = "R1"
      )
  }
  
  # read 2
  
  gcl_df <- NULL
  
  if (length(gcl) > 0) {
    
    gcl_df <- read_tsv(
      I(paste(gcl, collapse = "\n")),
      col_names = c(
        "type",
        "gc",
        "count"
      ),
      show_col_types = FALSE
    ) %>%
      mutate(
        sample = sample,
        read = "R2"
      )
  }
  
  bind_rows(
    gcf_df,
    gcl_df
  )
}

# Extract GC data
gc_data <- map_dfr(
  SAMPLES,
  extract_gc
)


# Normalize within each sample/read
gc_data <- gc_data %>%
  group_by(
    sample,
    read
  ) %>%
  mutate(
    proportion = count / sum(count)
  ) %>%
  ungroup()


# plot 

p_gc <- ggplot(
  gc_data,
  aes(
    x = gc,
    y = proportion,
    colour = sample
  )
) +
  geom_line(
    linewidth = 0.8
  ) +
  facet_wrap(
    ~ read,
    ncol = 1
  ) +
  scale_x_continuous(
    breaks = seq(0, 100, by = 10),
    minor_breaks = seq(0, 100, by = 1)
  ) +
  scale_y_continuous(
    labels = scales::percent,
    breaks = seq(0, 0.1, by = 0.02),
    minor_breaks = seq(0, 0.1, by = 0.01)
  ) +
  labs(
    title = paste(
      "GC-content distribution -",
      METHOD
    ),
    subtitle = "Normalized within each sample",
    x = "GC content (%)",
    y = "Proportion of reads",
    colour = "Sample"
  ) +
  theme_bw(
    base_size = 13
  ) +
  theme(
    plot.title = element_text(
      face = "bold"
    ),
    axis.title = element_text(
      face = "bold"
    ),
    panel.grid.minor = element_blank()
  )

print(p_gc)


ggsave(
  file.path(
    OUTDIR,
    paste0(
      "GC_distribution_R1_R2_",
      METHOD,
      ".png"
    )
  ),
  p_gc,
  width = 9,
  height = 8,
  dpi = 300
)


############################################
## mapping quality

extract_MQ <- function(sample) {
  
  lines <- read_stats(sample)
  
  # MAPQ section
  mq <- lines[
    grepl("^MAPQ\\t", lines)
  ]
  
  if (length(mq) == 0) {
    
    warning(
      "No MAPQ section found for sample ",
      sample
    )
    
    return(NULL)
  }
  
  mq_df <- read_tsv(
    I(paste(mq, collapse = "\n")),
    col_names = c(
      "type",
      "mapq",
      "count"
    ),
    show_col_types = FALSE
  ) %>%
    mutate(
      sample = sample
    )
  
  mq_df
}


# Extract MAPQ data
mq_data <- map_dfr(
  SAMPLES,
  extract_MQ
)


# Normalize within each sample
mq_data <- mq_data %>%
  group_by(
    sample
  ) %>%
  mutate(
    proportion = count / sum(count)
  ) %>%
  ungroup()


#plot

p_mq <- ggplot(
  mq_data,
  aes(
    x = mapq,
    y = proportion,
    colour = sample,
    linetype = sample
  )
) +
  geom_line(
    linewidth = 0.9
  ) +
  scale_y_continuous(
    labels = scales::percent
  ) +
  labs(
    title = paste(
      "Mapping quality distribution -",
      METHOD
    ),
    subtitle = "Normalized within each sample",
    x = "Mapping quality (MAPQ)",
    y = "Proportion of reads",
    colour = "Sample",
    linetype = "Sample"
  ) +
  theme_bw(
    base_size = 13
  ) +
  theme(
    plot.title = element_text(
      face = "bold"
    ),
    axis.title = element_text(
      face = "bold"
    ),
    panel.grid.minor = element_blank()
  )

print(p_mq)


ggsave(
  file.path(
    OUTDIR,
    paste0(
      "mapping_quality_distribution_",
      METHOD,
      ".png"
    )
  ),
  p_mq,
  width = 9,
  height = 6,
  dpi = 300
)
