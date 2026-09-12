library(dplyr)
library(readr)

# ------------------------------------------------------------
# Paths
# ------------------------------------------------------------

base_dir <- "/scratch/project_2019524/pinksalmon_alignments"

method_dirs <- c(
  lowinput = file.path(
    base_dir,
    "lowinput/coverage_depth/depth_stats"
  ),
  lowinput_P100 = file.path(
    base_dir,
    "lowinput_P100/coverage_depth/depth_stats"
  ),
  normalinput = file.path(
    base_dir,
    "normalinput/coverage_depth/depth_stats"
  )
)

# Output directory
out_dir <- file.path(
  base_dir,
  "depth_percent_comparison"
)

# Samples to include
samples <- c("8034", "8035", "8036", "8038")


# ------------------------------------------------------------
# Read all depth statistics files from all methods
# ------------------------------------------------------------

all_depth_stats <- lapply(names(method_dirs), function(method) {
  
  method_dir <- method_dirs[[method]]
  
  # Find all individual depth statistics files
  files <- list.files(
    method_dir,
    pattern = "\\.depth_stats\\.txt$",
    full.names = TRUE
  )
  
  lapply(files, function(file) {
    
    # Get sample ID from filename
    sample <- basename(file)
    sample <- sub("\\.depth_stats\\.txt$", "", sample)
    
    # Rename old sample IDs
    sample <- case_when(
      sample == "8167" ~ "8036",
      sample == "8168" ~ "8038",
      TRUE ~ sample
    )
    
    # Read file
    read.delim(
      file,
      header = TRUE,
      sep = "\t",
      stringsAsFactors = FALSE
    ) %>%
      mutate(
        Sample = sample,
        Method = method
      )
    
  })
  
}) %>%
  unlist(recursive = FALSE) %>%
  bind_rows()


# ------------------------------------------------------------
# Keep only desired samples
# ------------------------------------------------------------

all_depth_stats <- all_depth_stats %>%
  filter(Sample %in% samples)


# ------------------------------------------------------------
# Keep the required columns
# ------------------------------------------------------------

all_depth_stats <- all_depth_stats %>%
  select(
    Sample,
    Method,
    Total_bases,
    Mean_depth,
    Percent_5x,
    Percent_8x,
    Percent_10x,
    Percent_12x,
    Percent_15x
  )



for (sample in samples) {
  
  # Select this sample
  sample_data <- all_depth_stats %>%
    filter(Sample == sample)
  
  # Create sample directory
  sample_dir <- file.path(
    out_dir,
    sample
  )
  
  dir.create(
    sample_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  # Output filename
  output_file <- file.path(
    sample_dir,
    paste0("depth_stats_", sample, ".csv")
  )
  
  # Write CSV
  write_csv(
    sample_data,
    output_file
  )
  
  message("Written: ", output_file)
}
